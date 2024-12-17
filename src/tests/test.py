import psutil
import time
import sys
import os

# Adding root dir. of the project (OWLGENOME) to sys.path so Python can find the src module
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

# Import the required modules
#from src.index.index import minimize_hash, make_hash_table  # Updated after renaming import.py to inport
from src.alignment.nw_align import nw_align  # Import the actual function from nw_align.py
from src.alignment.sw_align import sw_align  # Import the actual function from sw_align.py

# Function to measure RAM usage, memory usage, iterations, and reads
def monitor_function(func, log=False, *args, **kwargs, ):
    '''
    Wrapper function to measure memory and time complexity of other functions
    Inputs: 
        func: function to monitor
        log: bool to determine whether values are printed
        *args: arguments of the monitored function
        **kwards: kwargs of the monitored function
    
    Output: 
        function output
        
    '''
    # Get the initial memory usage
    process = psutil.Process()
    mem_before = process.memory_info().rss / (1024 ** 2)  # in MB
    start_time = time.time()

    # Run function
    result = func(*args, **kwargs)

    ## Get end states
    end_time = time.time()
    mem_after = process.memory_info().rss / (1024 ** 2)  # in MB

    # Calculate performance metrics
    elapsed_time = end_time - start_time
    mem_usage = mem_after - mem_before

    if log:
        print(f"\nFunction: {func.__name__}")
        print(f"Memory before: {mem_before:.2f} MB, after: {mem_after:.2f} MB, used: {mem_usage:.2f} MB")
        print(f"Execution time: {elapsed_time:.4f} seconds")

    return result

def write_simulated_seqs(num_iter, genome_seq, seq_outpath, truth_outpath,
                        min_read_length = 500,
                        max_read_length = 12000,
                        revcom_prob = .5,
                        corruption_prob = 1,
                        error_rate = 0.15,
                        indel_rate = 0.05,
                        unmapped_seq_rate = 0.01,
                        ):
    '''
    Simulate sequences and write to fasta file
    '''
    ## Get indices for unmapped seqs
    unmapped_seq_idx = random.sample(range(num_iter), int(num_iter*unmapped_seq_rate))
    records = []
    for i in tqdm(range(num_iter)):
        if i not in unmapped_seq_idx:
            ## Extract query seq
            query_len = int((random.random() * (max_read_length - min_read_length)) + min_read_length)
            seq_start = int(random.random() * (len(genome_seq)-(query_len + 1)))
            seq = genome_seq[seq_start:seq_start+query_len]
            # Modify query seq
            if random.random() < revcom_prob:
                seq = RevCom(seq)
            if random.random() < corruption_prob:
                seq = corrupt_seq(seq, error_rate=error_rate, indel_rate=indel_rate)
            description = str(seq_start) + "-" + str(seq_start+query_len)
        else:
            query_len = int((random.random() * (max_read_length - min_read_length)) + min_read_length)
            seq = generate_dna(query_len)
            description = 'nan-nan'
        ## Convert to biopython SeqRecord
        record = SeqRecord(Seq(seq), id=str(i), description=description)
        records.append(record)
    
    ## Write sequences to file
    SeqIO.write(records, seq_outpath, 'fasta')

    ## Write ground truth
    locs = []
    for record in records: 
        id = record.id
        _start, end = record.description.split('-')
        _start = _start.split(' ')[-1]
        locs.append([id, _start, end])

    ## Write ground truth to file
    df = pd.DataFrame(locs)
    df.set_index(0, inplace=True)
    df.to_csv(truth_outpath, sep='\t', header=None)
    return 

def process_simulated_seqs(seq_path, seq_format):
    '''
    Reads in all simulated seqs and extracts posiiton information from descriptions
    '''
    records = [i for i in SeqIO.parse(seq_path, seq_format)]
    locs = []
    for record in records: 
        id = record.id
        _start, end = record.description.split('-')
        _start = _start.split(' ')[-1]
        locs.append((id, float(_start), float(end)))
    return locs

