import sys
import subprocess
import importlib.util

def check_package(name, import_name=None):
    if import_name is None:
        import_name = name
    
    try:
        spec = importlib.util.find_spec(import_name)
        if spec is None:
            return None
        
        module = importlib.import_module(import_name)
        return getattr(module, "__version__", "Installed (Unknown Version)")
    except ImportError:
        return None
    except Exception as e:
        return f"Error loading: {e}"

def get_cuda_version():
    try:
        output = subprocess.check_output(["nvcc", "--version"]).decode("utf-8")
        for line in output.split('\n'):
            if "release" in line:
                return line.split("release")[-1].split(",")[0].strip()
    except:
        return None

def check_requirements():
    print("="*60)
    print(f"Environment Verification: {sys.prefix}")
    print("="*60)

    # 1. Core System
    python_version = sys.version.split()[0]
    print(f"Python: {python_version}")
    
    cuda_sys = get_cuda_version()
    print(f"System CUDA (nvcc): {cuda_sys if cuda_sys else 'Not found'}")

    # 2. PyTorch & CUDA
    torch_ver = check_package("torch")
    print(f"\n[PyTorch]")
    print(f"  - torch: {torch_ver}")
    
    if torch_ver and torch_ver != "Error loading: No module named 'torch'":
        try:
            import torch
            print(f"  - CUDA Available: {torch.cuda.is_available()}")
            if torch.cuda.is_available():
                print(f"  - CUDA Device: {torch.cuda.get_device_name(0)}")
                print(f"  - Torch CUDA Version: {torch.version.cuda}")
                print(f"  - CUDA Architecture (sm): {torch.cuda.get_device_capability()}")
            else:
                 print("  - WARNING: PyTorch cannot see CUDA devices.")
        except Exception as e:
            print(f"  - Error checking CUDA: {e}")
    else:
        print("  - PyTorch not found or failed to load.")

    # 3. Common Packages
    print(f"\n[Common Packages]")
    common = ["numpy", "scipy", "pandas", "biopython", "matplotlib", "networkx", "requests", "tqdm"]
    for pkg in common:
        ver = check_package(pkg)
        print(f"  - {pkg}: {ver if ver else 'MISSING'}")

    # 4. Tool-Specific Checks
    print(f"\n[Tool Compatibility Checks]")
    
    # RFdiffusion 1 (SE3nv) requirements
    hydra = check_package("hydra", "hydra")
    dgl = check_package("dgl")
    pyrsistent = check_package("pyrsistent")
    
    rfd1_ready = all([hydra, dgl, pyrsistent, torch_ver])
    print(f"  RFdiffusion 1: {'[OK]' if rfd1_ready else '[MISSING]'} (Needs: hydra, dgl, pyrsistent)")

    # RFdiffusion 2 (protgen/cuda124) requirements
    e3nn = check_package("e3nn")
    pyrosetta = check_package("pyrosetta")
    openbabel = check_package("openbabel", "openbabel") # conda install often provides openbabel bindings
    
    rfd2_ready = all([hydra, dgl, e3nn, pyrosetta, torch_ver])
    print(f"  RFdiffusion 2: {'[OK]' if rfd2_ready else '[MISSING]'} (Needs: hydra, dgl, e3nn, pyrosetta)")

    # ProteinMPNN requirements
    # Minimal: torch, numpy
    pmpnn_ready = all([check_package("numpy"), torch_ver])
    print(f"  ProteinMPNN:   {'[OK]' if pmpnn_ready else '[MISSING]'} (Needs: numpy, torch)")

    # LigandMPNN requirements
    ml_coll = check_package("ml_collections")
    dm_tree = check_package("dm-tree", "tree")
    if not dm_tree: dm_tree = check_package("dm-tree", "dm_tree")
    prody = check_package("prody", "ProDy")
    if not prody: prody = check_package("prody", "prody")

    lmpnn_ready = all([ml_coll, dm_tree, prody, torch_ver])
    print(f"  LigandMPNN:    {'[OK]' if lmpnn_ready else '[MISSING]'} (Needs: ml_collections, dm-tree, prody)")

    # AllMetal3D requirements
    moleculekit = check_package("moleculekit")
    am3d_ready = all([moleculekit, torch_ver]) # Simplified check
    print(f"  AllMetal3D:    {'[OK]' if am3d_ready else '[MISSING]'} (Needs: moleculekit)")

    # SuperMetal requirements
    sklearn = check_package("scikit-learn", "sklearn")
    super_ready = all([sklearn, torch_ver])
    print(f"  SuperMetal:    {'[OK]' if super_ready else '[MISSING]'} (Needs: scikit-learn)")

    # 5. Graph Neural Networks (Common for RFd/LMPNN)
    print(f"\n[Graph Dependencies (PyG/DGL)]")
    print(f"  - dgl: {dgl if dgl else 'MISSING'}")
    
    pyg_pkgs = ["pyg_lib", "torch_scatter", "torch_sparse", "torch_cluster", "torch_spline_conv"]
    for pkg in pyg_pkgs:
        ver = check_package(pkg)
        print(f"  - {pkg}: {ver if ver else 'MISSING'}")

if __name__ == "__main__":
    check_requirements()
