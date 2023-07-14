import tempfile
import subprocess
import time

# collections.(Mapping,Sequence) were
# deprecated in newer version of python
try:
    from collections.abc import Sequence
except:
    from collections import Sequence

def multi_run(jobs, max_queue=None, verbose=False, submit_delay=2, update_delay=20):
  """Submit and monitor multiple SLURM jobs"""

  total = len(jobs)
  finished = 0
  running = 0
  active_jobs = []

  if max_queue is None:
    max_queue = total

  # Submit jobs
  while finished + running < total:
    while running < max_queue and finished + running < total:
      jobs[finished+running].submit()
      
      time.sleep(submit_delay)

      active_jobs.append(jobs[finished+running])
      running += 1

    time.sleep(update_delay)

    multi_update_status(active_jobs)

    active_jobs_updated = []
    for i in range(len(active_jobs)):
      if active_jobs[i].status in ["PENDING", "RUNNING"]:
        active_jobs_updated.append(active_jobs[i])
      else:
        #job errored out
        running -= 1
        finished += 1

    active_jobs = active_jobs_updated
  
  # Wait for jobs to finish
  while active_jobs:
    time.sleep(update_delay)

    multi_update_status(active_jobs)

    active_jobs_updated = []
    for i in range(len(active_jobs)):
      if active_jobs[i].status in ["PENDING", "RUNNING"]:
        active_jobs_updated.append(active_jobs[i])
      else:
        #job errored out
        running -= 1
        finished += 1

    active_jobs = active_jobs_updated

def multi_update_status(jobs, max_attempts=3, sleep_attempt=5):
  """"""
  for j in jobs:
    j.status = None

  attempt = 0
  while attempt < max_attempts and [j for j in jobs if j.status == None]:
    if attempt > 0:
      time.sleep(sleep_attempt)

    sacct_str = subprocess.check_output("sacct", shell=True)
    sacct_str = sacct_str.decode("utf-8")

    sacct_lines = sacct_str.split("\n")
    for line in sacct_lines[2:]:
      a = line.split()
      try:
        line_id = int(a[0])
      except:
        line_id = None
      
      for j in jobs:
        if line_id == j.id:
          j.status = a[4]
    
    attempt += 1

class job:
  def __init__(self, cmd, name, out_file, err_file, sb_file,
                queue='queue0', cpu=1, mem=None, time=None, 
                gpu=None, conda=None, modules=None):
    """
    """
    self.cmd = cmd
    self.name = name
    self.out_file = out_file
    self.err_file = err_file
    self.sb_file = sb_file
    self.queue = queue
    self.cpu = cpu
    self.mem = mem
    self.time = time
    self.gpu = gpu
    self.conda = conda

    if isinstance(modules, str):
      self.modules = [modules]
    elif isinstance(modules, Sequence):
      self.modules = modules
    else:
      self.modules = []

    self.id =  None
    self.status = None

  def submit(self):
    """
    run command on cluster
    """
    if self.sb_file is None:
      sb_tempfile = tempfile.NamedTemporaryFile()
      sb_file = sb_tempfile.name
    else:
      sb_file = self.sb_file
    
    with open(sb_file, "w") as sb_filehandle:

      sb_filehandle.write("#!/bin/bash\n\n")

      sb_filehandle.write(f"#SBATCH -p {self.queue}\n")
      sb_filehandle.write(f"#SBATCH -n 1\n")
      sb_filehandle.write(f"#SBATCH -c {self.cpu}\n")

      if self.name:
        sb_filehandle.write(f"#SBATCH -J {self.name}\n")
      if self.out_file:
        sb_filehandle.write(f"#SBATCH -o {self.out_file}\n")
      if self.err_file:
        sb_filehandle.write(f"#SBATCH -e {self.err_file}\n")
      if self.mem:
        sb_filehandle.write(f"#SBATCH --mem {self.mem}\n")
      if self.time:
        sb_filehandle.write(f"#SBATCH --time {self.time}\n")

      for module in self.modules:
        sb_filehandle.write(f"module load {module}\n")

      if self.conda:
        sb_filehandle.write(f"conda activate {self.conda}\n")

      sb_filehandle.write(self.cmd)
      sb_filehandle.write("\n")
    
    submit_str = subprocess.check_output(f"sbatch {sb_file}", shell=True)
    
    self.id = int(submit_str.split()[3])
  
  def update_status(self, max_attempts=3, sleep_attempt=5):
    """
    """
    status = None
    attempt = 0
    
    while attempt < max_attempts and status == None:
      if attempt > 0:
        time.sleep(sleep_attempt)
      sacct_str = subprocess.check_output('sacct', shell=True)
      sacct_str = sacct_str.decode('utf-8')

      sacct_lines = sacct_str.split("\n")
      for line in sacct_lines[2:]:
        a = line.split()
        try:
          line_id = int(a[0])
        except:
          line_id = None
        if line_id == self.id:
          status = a[5]

        attempt += 1

    if status == None:
      return False
    else:
      self.status = status
      return True
  

        
    

