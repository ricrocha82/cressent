#!/usr/bin/env python3

import os
import sys
import argparse
import importlib.util
import logging
from pathlib import Path

# Get the absolute path to the directory containing this script
SCRIPT_DIR = Path(__file__).resolve().parent
MODULES_DIR = SCRIPT_DIR / "modules"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("cressent")

def import_module_from_file(module_name, file_path):
    """Dynamically import a module from a file path."""
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    if spec is None:
        return None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module

def list_available_modules():
    """List all available modules in the modules directory."""
    modules = []
    for file_path in MODULES_DIR.glob("*.py"):
        if file_path.name.startswith("__"):
            continue
        module_name = file_path.stem
        modules.append(module_name)
    return sorted(modules)

def main():
    # Create the main parser
    parser = argparse.ArgumentParser(
        description="cressent: A comprehensive toolkit for ssDNA virus analysis",
        usage="cressent <module> [options]"
    )
    
    # Add the module argument
    parser.add_argument("module", help="Module to run")
    
    # If no arguments provided or only --help, list available modules
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] in ["-h", "--help"]):
        print("Available modules:")
        for module in list_available_modules():
            print(f"  {module}")
        parser.print_help()
        return 0
    
    # Parse just the module name
    args, remaining_args = parser.parse_known_args()
    module_name = args.module
    
    # Check if the module exists
    module_path = MODULES_DIR / f"{module_name}.py"
    if not module_path.exists():
        logger.error(f"Module '{module_name}' not found.")
        print("Available modules:")
        for module in list_available_modules():
            print(f"  {module}")
        return 1
    
    # Import the module
    try:
        module = import_module_from_file(module_name, module_path)
        
        # If the module has a 'main' function, call it with remaining arguments
        if hasattr(module, "main"):
            # Reset sys.argv for the module's argparse
            sys.argv = [sys.argv[0] + " " + module_name] + remaining_args
            return module.main()
        else:
            logger.error(f"Module '{module_name}' does not have a main function.")
            return 1
    except Exception as e:
        logger.error(f"Error running module '{module_name}': {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())