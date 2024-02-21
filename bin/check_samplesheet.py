#!/usr/bin/env python

"""Provide a command line tool to validate the structure of tabular samplesheets."""

import argparse
import csv
import logging
import sys
from pathlib import Path

logger = logging.getLogger()


REQUIRED_COLUMNS = frozenset({"sample", "data_directory", "n_cell_types", "bleeding_correction", "spatial_transcriptional_programs"})

class RowChecker:
    """
    Define a service that can validate the structure of each given row.

    Attributes:
        rows (list): A list of dicts, each dict corresponds to a row in the file.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.rows = []

    def check_row_structure(self, row):
        """
        Check if the row has the required structure.
        In this case, we are only checking if the two required columns are present.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).
        """
        self.rows.append(row)


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the required columns.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet.
        file_out (pathlib.Path): Where the validated samplesheet should be created.

    """

    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle)
        # Validate the existence of the expected header columns.
        if not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            req_cols = ", ".join(REQUIRED_COLUMNS)
            logger.critical(f"The sample sheet must contain these column headers: {req_cols}.")
            sys.exit(1)

        # Validate each row.
        checker = RowChecker()
        for row in reader:
            checker.check_row_structure(row)

    # Writing rows to the output file.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, REQUIRED_COLUMNS, delimiter=",")
        writer.writeheader()
        writer.writerows(checker.rows)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate the structure of a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv validated_samplesheet.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Input samplesheet in CSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Output validated samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
