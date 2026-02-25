#!/usr/bin/env python3
import sys
import subprocess
import os

def check_modules():
    print("Testing Foundry Environment Setup...")
    
    try:
        import torch
        print(f"✅ PyTorch is working. Version: {torch.__version__}")
        print(f"✅ CUDA available: {torch.cuda.is_available()}")
    except ImportError:
        print("❌ PyTorch is not installed or failed to import.")
        sys.exit(1)

    print("\nTesting RosettaCommons Foundry Imports...")
    modules_to_test = [
        ("foundry", "Foundry core package"),
        ("rfd3", "RFdiffusion(3)"),
        ("mpnn", "LigandMPNN"),
        ("rf3", "RoseTTAFold(3)"),
    ]

    all_passed = True
    for mod_name, desc in modules_to_test:
        try:
            # Check if module can be imported
            __import__(mod_name)
            print(f"✅ Successfully imported {desc} ({mod_name})")
        except ImportError as e:
            print(f"❌ Failed to import {desc} ({mod_name}): {e}")
            all_passed = False

    if all_passed:
        print("\n✨ All components successfully verified!")
    else:
        print("\n⚠️  Some modules failed to import. Please check your installation.")
        sys.exit(1)

def install_models():
    print("\nChecking for Foundry Models...")
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    print(f"Using checkpoint directory: {checkpoint_dir}")

    try:
        # Check what models are currently installed
        result = subprocess.run(
            ["foundry", "list-installed"],
            capture_output=True,
            text=True,
            check=True
        )
        installed_output = result.stdout.lower()
        
        models_to_check = ["rfd3", "rf3", "ligandmpnn"]
        all_installed = all(m in installed_output for m in models_to_check)
        
        if all_installed:
            print("✅ Base models are already installed.")
        else:
            print("📥 Base models not fully found. Installing...")
            try:
                subprocess.run(["foundry", "install", "base-models", "--checkpoint-dir", checkpoint_dir], check=True)
                print("✅ Successfully installed base models.")
            except subprocess.CalledProcessError:
                print("❌ Failed to download base models.")
                    
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to run foundry CLI to check installed models. Error: {e}")
    except FileNotFoundError:
        print("❌ Foundry CLI not found. Is the environment activated?")

if __name__ == "__main__":
    check_modules()
    install_models()
