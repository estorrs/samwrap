import os
import subprocess
import sys
from collections import OrderedDict
from multiprocessing import Pool

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_sort_command(bam_fp, output_fp, sort_threads=1):
    """get samtools sort command

    arguments

    bam_fp
        - input bam filepath. set to None if input is to be streamed.
    output_fp
        - filepath for sorted output bam. set to None to stream output to std out.
    """
    samtools_command = ['samtools', 'sort']

    if output_file is not None:
        samtools_command += ['-o', output_fp]

    if bam_fp is not None:
        samtools_command.append(bam_fp)

    return samtools_command

def get_position_filter_command(bam_fp, output_fp, positions_fp):
    """get samtools position filter command

    arguments

    bam_fp
        - input bam filepath. set to None if input is being streamed
    output_fp
        - filepath for sorted output bam. set to None to stream output to std out.
    positions_fp
        - filepath for positions .bed file. File contains tab seperated positions in the following format:
        <chrom>\t<start-pos>\t<stop-pos>
    """
    samtools_command = ['samtools', 'view', '-L', positions_fp]

    if output_file is not None:
        samtools_command += ['-o', output_fp]

    if bam_fp is not None:
        samtools_command.append(bam_fp)

    return samtools_command

def index_bam(bam_fp):
    """Index the given bam file.

    The indexed .bai file will be put in same directory as input file
    """
    samtools_command = ['samtools', 'index', bam_fp]
    subprocess.check_output(samtools_command)


class SamtoolsWrapper(object):
    def __init__(self, input_files, output_dir, operations_dict,
                 output_descriptor='samwrap', threads=1, verbose=True):
        """Wrapper for Samtools
        
        args:
        
        input_files: list
            list of input bams
        output_dir: str
            directory to store outputs in
        operations_dict: OrderedDict
            dict storing the samtools operations to perform on the input files.
            operations will be preformed in the same order as in the ordered dict.
            
            {<operation_identifier>: <operation_kwargs>}
            
            possible identifiers - "index", "sort", "position_filter"
            
            operation_kwargs - arguments for the operation.
                for position filter - {'positions_fp': '/path/to/positions.bed'}
                for sort - {'sort_threads': num_sort_threads}
                for index - None
            
        kwargs:
        
        threads: int
            number of precesses to use
        verbose: bool
            verbose output
            
        """
        self.input_files = input_files
        self.output_dir = output_dir
        
        self.operations_dict = operations_dict
        
        self.output_descriptor = output_descriptor
        self.num_threads = threads
        self.verbose = verbose
        
    def run_bams(self):
        # pull out index if it's in operations
        v = self.operations_dict.get('index')
        if v is not None:
            index_input_bam = True
            del self.operations_dict['index']
        else:
            index_input_bam = False
        
        # create argument pool
        arg_pool = []
        for fp in self.input_files:
            sample = fp.split('/')[-1]
            sample = re.sub(r'.bam', self.output_descriptor + '.bam', sample)
            output_fp = os.path.join(self.output_dir, sample)
            arg_pool.append((fp, output_fp, index_input_bam))
        
        with Pool(self.threads) as p:
            p.map(self.worker_wrapper, arg_pool)
        
    def worker_wrapper(self, args):
        return self.execute_bam(*args)
        
    def execute_bam(self, bam_fp, output_fp, index_input_bam=False):
        """Execute a bam file"""
        # index bam if needed
        if index_input_bam:
            index_bam(bam_fp)
        
        # execute remaining commands
        command_args = self.generate_stream_commands(bam_fp, output_fp)
        subprocess.check_output(command_args)
        
        if self.verbose:
            eprint(f'{bam_fp} completed')
        
        return output_fp
    
    def generate_stream_commands(self, bam_fp, output_fp):
        """
        return command to run input bam through the given operations
        
        """
        command_args = []
        
        items = self.operations_dict.items()
        
        # get first command
        operation_identifier, operation_kwargs = next(items)
        num_operations = len(self.operations_dict)
        # if only one operation return output
        if num_operations == 0:
            return self.get_command_from_identifier(bam_fp, output_fp, operation_identifier, operation_kwargs)
        else:
            command_args += self.get_command_from_identifier(bam_fp, None, operation_identifier, operation_kwargs) + ['|']
            
        c = 0
        for operation_identifier, operation_kwargs in items:
            # if last operation, write an output file, otherwise stream to next operation
            if num_operations - 1 == c:
                command_args += self.get_command_from_identifier(None, output_fp, operation_identifier, operation_kwargs)
            else:
                command_args += self.get_command_from_identifier(None, None, operation_identifier, operation_kwargs) + ['|']
        
        return command_args
    
    @staticmethod
    def get_command_from_identifier(bam_fp, output_fp, identifier, kwargs):
        """return command given identifier and operation kwargs"""
        if identifier == 'sort':
            return get_sort_command(bam_fp, output_fp, sort_threads=kwargs['sort_threads'])
        if identifier == 'position_filter':
            return get_position_filter_command(bam_fp, output_fp, kwargs['positions_fp'])

