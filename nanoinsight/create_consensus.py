import os
#import psutil
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
# from Bio.Align.Applications import MafftCommandline
from Bio.motifs import Motif
from concurrent.futures import ThreadPoolExecutor
import subprocess
import shutil
################################################################################################################################################################
#(2) Generating multiple alignment and consensus sequences for each INS
################################################################################################################################################################              

def create_cons(vcf, wk_dir, fasta_dir, id_seq, num_threads, mafft_path, batch_size=100, num_parallel_workers=5):
    """
    # # Get the total number of CPU threads  ##### How about just allow the user to specify number of threads to use?
    # total_cpus = os.cpu_count()
    # # Calculate the available CPU threads based on CPU usage
    # cpu_usage = psutil.cpu_percent(interval=1)  
    # available_cpus = total_cpus * (1 - (cpu_usage / 100))
    # # Calculate num_threads as 75% of the available CPUs
    # num_threads = int(available_cpus * 0.75)
    # # Ensure the number of threads is at least 1
    # num_threads = max(num_threads, 1)
    #print ('No of threads is used for Mafft is', num_threads)
    # Define your batch size and number of parallel workers
    # batch_size = 100
    # num_parallel_workers = 5  #process 5 batches concurrently
    """
    threads_per_job=int(num_threads/num_parallel_workers)  ##### What if num_threads is < 5? 
    #print ('No of threads for each job is', threads_per_job)
    #fasta_dir = "./fasta_files"
    MA_dir = os.path.join(wk_dir, "MA")
    con_dir = os.path.join(wk_dir, "con")
    mafft_exe = mafft_path #input arg
    os.makedirs(MA_dir, exist_ok=True)
    os.makedirs(con_dir, exist_ok=True)
    threads_per_job = submit_jobs(fasta_dir, MA_dir, con_dir, mafft_exe, threads_per_job, num_parallel_workers, batch_size)
    rename_header(con_dir, id_seq)
    con_fasta = cat_consensus(vcf, con_dir, wk_dir)
    remove_tmpdirs(fasta_dir, MA_dir, con_dir)
    return con_fasta, threads_per_job

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
                # mafft_cline = MafftCommandline(mafft_exe,input=input_path, auto=True, thread=threads_per_job)
                # stdout, stderr = mafft_cline()
                # with open(MA_file, "w") as out_file:
                #     out_file.write(stdout)
                mafftcmdline(mafft_exe, input_path, MA_file, threads_per_job)
                # Generate consensus sequence for the output file
                # alignment = AlignIO.read(MA_file, "fasta")
                # summary = AlignInfo.SummaryInfo(alignment)
                # consensus = summary.dumb_consensus(threshold=0.51, ambiguous='X')
                msa = AlignIO.read(MA_file, "fasta")
                alignment = msa.alignment
                motif = Motif('atgc', alignment)
                counts = motif.counts
                consensus = counts.calculate_consensus(identity=0.51)
                # Write consensus sequence to output file
                consensus_record = SeqRecord(Seq(consensus), id="consensus")
                SeqIO.write(consensus_record, con_file, "fasta")
    
    files = [f for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
    with ThreadPoolExecutor(max_workers=num_parallel_workers) as executor:
        for i in range(0, len(files), batch_size):
            batch = files[i:i + batch_size]
            executor.submit(process_batch, batch)
    return threads_per_job

## Header of each file should be renamed to the SV ID
def rename_header(con_dir, id_seq):
    for con_file in os.listdir(con_dir):
        if con_file.endswith(".fasta"):
            # Split filename and extension
            name, ext = os.path.splitext(con_file)
            # Remove the '.con' extension from the name
            name = name.replace('.con', '')
            # Add coordinates of SV to the header line
            # matching_row = id_seq[id_seq['ID'] == name]
            # SV_info = matching_row['sv_coo'].values[0]
            # Open file and read lines
            with open(os.path.join(con_dir, con_file), "r") as f:
                lines = f.readlines()
                new_header = f">{name}\n"
            # Replace first line with SV ID and coordinates info if or only SV ID if header length > 50 
                #header_info = f"{name}{SV_info}"
                # if len(header_info) > 47:
                #     new_header = f">{name}\n"
                # else:
                #     new_header = f">{header_info}\n"
                lines[0] = new_header
            # Open file for writing and write updated lines
            with open(os.path.join(con_dir, con_file), "w") as f:
                f.writelines(lines)

# Concatenate all consensus sequences in one file for RM
def cat_consensus(vcf, con_dir, wk_dir):
    filename, extension = os.path.splitext(os.path.basename(vcf))
    samplename = filename.replace('.nanovar.pass', '')
    if '.vcf' in samplename:
        samplename = samplename.replace('.vcf', '')
    #samplename = os.path.basename(vcf).strip('.nanovar.pass.vcf.gz')
    out = os.path.join(wk_dir, f'{samplename}.ins.con.fasta')
    subprocess.run([f"cat {con_dir}/*.fasta > {out}"], capture_output=True, text=True, shell = True)
    return out

def remove_tmpdirs(fasta_dir, MA_dir, con_dir):
    tmp_dirs = [fasta_dir, MA_dir, con_dir]
    for tmp_dir in tmp_dirs:
        shutil.rmtree(tmp_dir)

# Run MAFFT command line
def mafftcmdline(mafft_exe, input, output, threads):
    p = subprocess.run([f"{mafft_exe} --auto --thread {threads} {input} > {output}"],
                       capture_output=True, text=True, shell=True)
    if p.returncode != 0:
        err = p.stderr.strip()
        raise Exception('Mafft returned non-zero exit status:\n %s' % err)
