import os

from seq_deposit.logger import get_logger

log = get_logger("util")


def create_directories() -> None:
    """Create the necessary directories for the output.

    Returns:
        None
    """
    directories = [
        "seq-deposit-output/dna",
        "seq-deposit-output/rna",
        "seq-deposit-output/fastas",
    ]
    if not os.path.exists("seq-deposit-output"):
        os.makedirs("seq-deposit-output", exist_ok=True)
    log.info("Outputs will be written in seq-deposit-output")
    for directory in directories:
        os.makedirs(directory, exist_ok=True)
    log.info("Directories created: %s", directories)
