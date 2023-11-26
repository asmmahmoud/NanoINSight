import os
#import psutil
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from concurrent.futures import ThreadPoolExecutor
import subprocess

################################################################################################################################################################
#(2) Generating multiple alignment and consensus sequences for each INS
################################################################################################################################################################              

def create_cons(vcf, wk_dir, fasta_dir, id_seq, num_threads, mafft_path, batch_size=100, num_parallel_workers=5):
    """
    # # Get the total number of CPU threads  ##### How about just allow the user to specify number of threads to use?
    # total_cpus = os.cpu_count()
    # # Calculate the available CPU threads based on CPU usage
    # cpu_usage = psutil.cpu_percent(interval=1)  # Adjust the interval as needed
    # available_cpus = total_cpus * (1 - (cpu_usage / 100))
    # # Calculate num_threads as 75% of the available CPUs
    # num_threads = int(available_cpus * 0.75)
    # # Ensure the number of threads is at least 1
    # num_threads = max(num_threads, 1)
    #print ('No of threads is used for Mafft is', num_threads)
    # Define your batch size and number of parallel workers
    # batch_size = 100
    # num_parallel_workers = 5  # You want to process 5 batches concurrently
    """
    threads_per_job=int(num_threads/num_parallel_workers)  ##### What if num_threads is < 5? 
    #print ('No of threads for each job is', threads_per_job)
    #fasta_dir = "./fasta_files"
    MA_dir = os.path.join(wk_dir, "MA")
    con_dir = os.path.join(wk_dir, "con")
    mafft_exe = mafft_path #input arg
    os.makedirs(MA_dir, exist_ok=True)
    os.makedirs(con_dir, exist_ok=True)
    submit_jobs(fasta_dir, MA_dir, con_dir, mafft_exe, threads_per_job, num_parallel_workers, batch_size)
    rename_header(con_dir, id_seq)
    cat_consensus(vcf, con_dir, wk_dir)

## Run Mafft and generate consensus with multi-threading
def submit_jobs(fasta_dir, MA_dir, con_dir, mafft_exe, threads_per_job, num_parallel_workers, batch_size):
  
    def process_batch(batch):
        for input_file in batch:
            if input_file.endswith(".fasta"):
                # Input and output file paths
                input_path = os.path.join(fasta_dir, input_file)
                MA_file = os.path.join(MA_dir, input_file.replace(".fasta", ".MA.fasta"))
                con_file = os.path.join(con_dir, input_file.replace(".fasta", ".con.fasta"))
                # Run MAFFT with MafftCommandline
                mafft_cline = MafftCommandline(mafft_exe,input=input_path, auto=True, thread=threads_per_job)
                stdout, stderr = mafft_cline()
                with open(MA_file, "w") as out_file:
                    out_file.write(stdout)
                # Generate consensus sequence for the output file
                alignment = AlignIO.read(MA_file, "fasta")
                summary = AlignInfo.SummaryInfo(alignment)
                consensus = summary.dumb_consensus(threshold=0.51, ambiguous='X')
                # Write consensus sequence to output file
                consensus_record = SeqRecord(Seq(consensus), id="consensus")
                SeqIO.write(consensus_record, con_file, "fasta")
    
    files = [f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    with ThreadPoolExecutor(max_workers=num_parallel_workers) as executor:
        for i in range(0, len(files), batch_size):
            batch = files[i:i + batch_size]
            executor.submit(process_batch, batch)
    #print ('Generating consensus sequences is DONE')

## Header of each file should be renamed to the SV ID
def rename_header(con_dir, id_seq):
    for con_file in os.listdir(con_dir):
        if con_file.endswith(".fasta"):
            # Split filename and extension
            name, ext = os.path.splitext(con_file)
            # Remove the '.con' extension from the name
            name = name.replace('.con', '')
            # Add coordinates of SV to the header line
            matching_row = id_seq[id_seq['ID'] == name]
            info_to_add = matching_row['sv_coo'].values[0]
            # Open file and read lines
            with open(os.path.join(con_dir, con_file), "r") as f:
                lines = f.readlines()
            # Replace first line with filename (without extension) and remove ".con.fasta" if it's there
            new_header = f">{name}{info_to_add}\n"
            lines[0] = new_header
            # Open file for writing and write updated lines
            with open(os.path.join(con_dir, con_file), "w") as f:
                f.writelines(lines)
    #print ('Renaming is DONE')

# Concatenate all consensus sequences in one file for RM
def cat_consensus(vcf, con_dir, wk_dir):
    samplename, ext = os.path.splitext(vcf)
    samplename = samplename.replace('.nanovar.pass', '')
    out = os.path.join(wk_dir, f'{samplename}.ins.con.fasta')
    subprocess.run([f"cat {con_dir}/*.fasta > {out}"], capture_output=True, text=True, shell = True)