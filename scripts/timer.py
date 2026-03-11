import subprocess
import resource
import time
import sys

start_time = time.time()
subprocess.run(sys.argv[1:], shell=True)
end_time = time.time()

# Get peak memory usage in KB
usage = resource.getrusage(resource.RUSAGE_CHILDREN)
print(f"\n--- Metrics ---")
print(f"Runtime: {end_time - start_time:.2f} seconds")
print(f"Peak Memory: {usage.ru_maxrss} KB")