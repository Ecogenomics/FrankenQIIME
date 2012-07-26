import logging
import tempfile
import subprocess
import pysam
from qiime.assign_taxonomy import TaxonAssigner
import os
import sys

class BwaswTaxonAssigner(TaxonAssigner):
    NO_HIT_STRING = 'No BWA-SW hit'
    
    def __init__(self, params):
        """ Initialize the object
        """
        _params = {
            'Threads': 10,
            }
        _params.update(params)
        TaxonAssigner.__init__(self, _params)
    
    def __call__(self, seq_path=None, seqs=None, result_path=None, log_path=None):
        """Returns dict mapping {seq_id:(taxonomy, confidence)} for each seq.
        """
        assert seq_path, \
         "Must provide seq_path when calling a BwaTaxonAssigner."
         
        # initialize the logger
        logger = self._get_logger(log_path)
        logger.info(str(self))
        
        # assign the bwasw index
        db = self.Params['database']
        
        # build the mapping of sequence identifier 
        # (wrt to the blast db seqs) to taxonomy
        id_to_taxonomy_map = self._parse_id_to_taxonomy_file(\
         open(self.Params['id_to_taxonomy_filepath'],'U'))
        eg_key = id_to_taxonomy_map.keys()[0]
        if eg_key is None:
            logger.error("Unexpectedly found 0 entries in the sequence id to taxonomy map file '%s'",self.Params['id_to_taxonomy_filepath'])
        else:
            logger.debug("Read %s entries of an ID to taxonomy map e.g. '%s' => '%s'",eg_key,id_to_taxonomy_map[eg_key])
            
        # Assign taxonomy to the sequences and return
        best_hits = self._seqs_to_best_hits(\
         seq_path, db, logger=logger)
        logger.info("Parsed %s hits from BWASW output",len(best_hits))
        
        result = self._map_ids_to_taxonomy(best_hits, id_to_taxonomy_map)
        
        # Write log data if we have a path (while the logger can handle
        # being called if we are not logging, some of these steps are slow).
        if log_path is not None:
            num_inspected = len(result)
            logger.info('Number of sequences inspected: %s' % num_inspected)
            num_null_hits = [r for r in result.values()].count(None)
            logger.info('Number with no blast hits: %s' % num_null_hits)
        
        if result_path:
            # if the user provided a result_path, write the 
            # results to file
            of = open(result_path,'w')
            for seq_id, result_dict in result.items():
                if result_dict['TAXONOMY'] == self.NO_HIT_STRING:
                    of.write('%s\t%s\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\n' %
                    (seq_id, result_dict['TAXONOMY']))
                else:
                    of.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % 
                        (seq_id,
                        result_dict['TAXONOMY'],
                        '-',#blast_result_dict['E-VALUE'],
                        result_dict['SUBJECT ID'],
                        '-',#blast_result_dict['% IDENTITY'],
                        '-',#blast_result_dict['ALIGNMENT LENGTH'],
                        '-',#blast_result_dict['MISMATCHES'],
                        '-',#blast_result_dict['GAP OPENINGS'],
                        '-',#blast_result_dict['Q. START'],
                        '-',#blast_result_dict['Q. END'],
                        '-',#blast_result_dict['S. START'],
                        '-')#blast_result_dict['S. END'])
                        )
            of.close()
            result = None
            logger.info('Result path: %s' % result_path)
        else:
            # Returning the data as a dict, so no modification to result
            # is necessary.
            pass
                 
            # if no result_path was provided, return the data as a dict
            logger.info('Result path: None, returned as dict.')
        #
        ## clean-up temp blastdb files, if a temp blastdb was created
        #if 'reference_seqs_filepath' in self.Params:
        #    map(remove,db_files_to_remove)
        #
        
        return result
    
    def _seqs_to_best_hits(self, seq_path, database, logger=None):
        """ Assign taxonomy to (seq_id,seq) pairs
        Returns dict mapping {seq_id:(taxonomy, confidence)} for each seq,
        by running BWA-SW (querying seq_path against the database), and then parsing
        the output file and mapping them through the id_to_taxonomy_map,
        by calling _get_taxonomy_identifers_from_sam_file().
        """
        
        # Use BWASW
        # Create a temporary directory to store results in
        sam_file_tuple = tempfile.mkstemp()
        sam_file_fd = sam_file_tuple[0]
        sam_file_path = sam_file_tuple[1]
        
        # Create stderr file for bwa
        stderr_tuple = tempfile.mkstemp()
        stderr_fd = stderr_tuple[0]
        stderr_path = stderr_tuple[1]
        
        num_threads = str(self.Params['Threads'])
        
        worked_ok = False
        try:
            # bwa bwasw database.fasta long_read.fastq > aln.sam
            # TODO: when bwa is not found in the path, it errors out with a not very useful error
            bwa_arguments = ['bwa','bwasw', '-t', num_threads, database, seq_path]
            if logger is not None:
                logger.info("Running BWA-SW like this: %s", bwa_arguments)
            return_code = subprocess.check_call(\
             bwa_arguments,\
             stdout = sam_file_fd,\
             stderr = stderr_fd
            )
            if logger is not None:
                logger.debug("BWA-SW return code: %i",return_code)
            
            # Convert a samfile into a list of sequence identifiers
            tax_ids = self._get_taxonomy_identifers_from_sam_file(sam_file_path,logger=logger)
            worked_ok = True
            
        except subprocess.CalledProcessError as e:
            bwa_stderr = open(stderr_path).read()
            if logger is not None:
                logger.error("There was a problem running BWA-SW: %s",e)
                logger.error("The STDERR was: %s",bwa_stderr)
            error = e
        
        finally:
            # Always remove the tempfile
            os.remove(sam_file_path)
            os.remove(stderr_path)
        
        
        if worked_ok:
            return tax_ids
        
        else:
            if logger is not None:
                logger.error("Some error was detected running BWA, failing.")
            sys.stderr.write(bwa_stderr)
            raise error
    
    def _get_taxonomy_identifers_from_sam_file(self, sam_file_path, logger=None):
        """ Take a sam file, and return a dictionary of query identifiers to identifiers from
        the database used in the creation of the sam file.
        
        BWA-SW detects chimeras and reports this by including two separate
        entries in the output SAM file. Here this is handled by assigning
        the sequence to the first hit in the SAM file, and logging that a
        chimera was detected to the logger.
        """
        
        samfile = pysam.Samfile(sam_file_path,'r')

        tax_identifiers = {}
        # For each read, add the hit identifier as the taxon identifier to the returned array
        for read in samfile.fetch(until_eof=True):
            query_name = read.qname
            
            if read.is_unmapped:
                tax_identifiers[query_name] = {}
                tax_identifiers[query_name]['SUBJECT ID'] = self.NO_HIT_STRING
            
            else:
                ref_name = samfile.getrname(read.tid)
                # if there is already a taxonomy assigned for this tag, bwasw detected a chimera
                if query_name in tax_identifiers.keys():
                    tax_identifiers[query_name]['CHIMERIC'] = True
                    if logger is not None:
                        logger.warn("BWA-SW flagged a tag as chimeric: '%s'", query_name)
                else:
                    # regular everyday hit to a sequence in the database
                    tax_identifiers[query_name] = {}
                    tax_identifiers[query_name]['SUBJECT ID'] = ref_name
        
        return tax_identifiers


    def _get_logger(self, log_path=None):
        if log_path is not None:
            handler = logging.FileHandler(log_path, mode='w')
        else:
            class NullHandler(logging.Handler):
                def emit(self, record): pass
            handler = NullHandler()
        logger = logging.getLogger("BwaswTaxonAssigner logger")
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _map_ids_to_taxonomy(self, hits, id_to_taxonomy_map):
        """ map {query_id:{best_blast_info_dict}} to {query_id:{best_blast_info_and_tax_dict}}
        """
        for query_id, hit_info in hits.items():
            query_id=query_id.split()[0]
            try:
                hits[query_id]['TAXONOMY'] = id_to_taxonomy_map.get(hit_info['SUBJECT ID'], self.NO_HIT_STRING)
            except TypeError:
                hits[query_id] = {'TAXONOMY': self.NO_HIT_STRING}

        return hits

