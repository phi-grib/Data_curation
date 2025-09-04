import subprocess
import sys
import os
import shutil
import pandas as pd

# Define the CLI command
CLI_COMMAND = "datacur"
TEST_ENDPOINT = "test-endpoint"
TEST_FILE = "sample_smiles.xlsx"

def run_command(command, expected_exit_code=0):
    """
    Runs a shell command and checks its exit code.
    Returns the command's stdout.
    """
    print(f"\nRunning command: {' '.join(command)}")
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        print("Success!")
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}.")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise

def setup_test_environment():
    """
    Sets up a clean test environment by configuring the tool in silent mode.
    """
    print("--- Setting up test environment ---")
    try:
        # Use a silent configuration to avoid user prompts
        run_command([CLI_COMMAND, "-c", "config", "-a", "silent"])
    except Exception as e:
        print("Setup failed. Please ensure the tool is installed and 'curate' is a recognized package name in your pyproject.toml.")
        raise e

def cleanup_test_environment():
    """
    Removes the test endpoint and configuration to clean up after the tests.
    """
    print("\n--- Cleaning up test environment ---")
    try:
        # Use the remove action to delete the test endpoint
        run_command([CLI_COMMAND, "-c", "manage", "-e", TEST_ENDPOINT, "-a", "remove"])
        # Remove the generated output file
        if os.path.exists("curation.tgz"):
            os.remove("curation.tgz")
        if os.path.exists("datacur.log"):
            os.remove("datacur.log")
    except Exception as e:
        print("Cleanup failed. Manual cleanup may be required.")
        print(f"Error: {e}")

def main():
    """
    Main function to run all tests.
    """
    try:
        setup_test_environment()

        # Test 1: Check if the help command works
        print("\n--- Test 1: Checking help command ---")
        stdout = run_command([CLI_COMMAND, "-h"])
        if "usage: datacur" in stdout:
            print("Test 1 Passed: Help message found.")
        else:
            print("Test 1 Failed: Help message not found.")
            return

        # Test 2: Test creating a new endpoint
        print("\n--- Test 2: Creating a new endpoint ---")
        run_command([CLI_COMMAND, "-c", "manage", "-a", "new", "-e", TEST_ENDPOINT])
        print("Test 2 Passed: New endpoint created successfully.")

        # Test 3: Test listing endpoints
        print("\n--- Test 3: Listing endpoints ---")
        stdout = run_command([CLI_COMMAND, "-c", "manage", "-a", "list"])
        if TEST_ENDPOINT in stdout:
            print("Test 3 Passed: Endpoint found in list.")
        else:
            print("Test 3 Failed: Endpoint not found in list.")
            return

        # Test 4: Curate a sample Excel file
        print("\n--- Test 5: Curating a sample file ---")
        run_command([CLI_COMMAND, "-i", TEST_FILE, "-e", TEST_ENDPOINT, "-c", "curate", "-a", "chem", "-s", "structure", "-id", "name", "-r"])
        print("Test 5 Passed: Curation command ran successfully.")
        
        # Test 5: Download the curated file
        print("\n--- Test 6: Downloading the curated file ---")
        run_command([CLI_COMMAND, "-c", "manage", "-a", "download", "-e", TEST_ENDPOINT, "-f", "sdf"])
        
        # Check if the file was created
        if os.path.exists("curation.tgz"):
            print("Test 6 Passed: Downloaded file found.")
        else:
            print("Test 6 Failed: Downloaded file not found.")
            return

        print("\n--- All basic installation and functionality tests completed successfully! ---")

    except Exception:
        print("\n--- One or more tests failed. ---")
    finally:
        cleanup_test_environment()

if __name__ == "__main__":
    main()