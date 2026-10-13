"""psutil stand-in for the Pyodide build: pyhmmer only calls cpu_count, and WebAssembly runs single-threaded."""
def cpu_count(logical=True):
    return 1
