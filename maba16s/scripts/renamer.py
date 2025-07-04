"""
Rename barcode directories in `input_directory` according to mapping in `samplesheet`.

Parameters:
    input_directory (str): Path to directory containing barcode folders.
    spreadsheet (str): Path to CSV, TSV, or Excel file with 2 columns:
                       [sample_name, barcode_name].

The code checks whether the provided sample sheet is valid and has the option of sending emails
should the sample sheet not adhere to the specified requirements. If no email addresses are provided
as arguments this functionality is skipped.
"""

import os
import shutil
import pandas as pd
import re
import sys
import smtplib
from email.message import EmailMessage
import argparse

### Function to send error emails ###
def send_error_email(subject, body, recipients):
    if not recipients:
        return

    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = f"{os.getlogin()}@localhost"  
    msg['To'] = ', '.join(recipients)
    msg.set_content(body)

    try:
        with smtplib.SMTP('localhost') as s:
            s.send_message(msg)
    except Exception as e:
        print(f"ERROR: Failed to send email: {e}", file=sys.stderr)

### Function to load and chec sample sheet
def validate_samplesheet(samplesheet, technician, mmb):
    
    recipients = [email for email in [technician, mmb] if email]

    # Sample sheet file check basics
    try:
        df = pd.read_excel(samplesheet)
    except Exception as e:
        msg = f"ERROR: Unable to read sample sheet - please check the file'{samplesheet}': {e}"
        send_error_email("MABA16S error", msg, recipients)
        sys.exit(msg)

    if df.shape[1] < 2:
        msg = "ERROR: The provided sample sheet must have at least 2 columns"
        send_error_email("MABA16S error", msg, recipients)
        sys.exit(msg)

    # Check the column names
    df.columns = [col.strip().lower() for col in df.columns]

    expected_cols = ['sample', 'barcode']
    if df.columns[0] != expected_cols[0] or df.columns[1] != expected_cols[1]:
        msg = (f"ERROR: The first two columns in the sample sheet must be named 'sample' and 'barcode'.\n"
               f"Found: {df.columns[0]} and {df.columns[1]}.\n"
               f"Please check the sample sheet and try again")
        send_error_email("MABA16S error", msg, recipients)
        sys.exit(msg)
    
    # Check the column content
    df = df[['sample', 'barcode']].dropna()

    # Ensure that sample names are unique
    if not df['sample'].is_unique:
        dupes = df['sample'][df['sample'].duplicated()].tolist()
        msg = (f"ERROR: Sample names in the sample sheet must be unique.\n"
               f"These samples are duplicated: {dupes}.\n"
               f"Please change the sample names and try again")
        send_error_email("MABA16S error", msg, recipients)
        sys.exit(msg)

    # Check that the barcodes are specified correctly
    df['barcode'] = df['barcode'].astype(str).str.strip().str.lower()

    invalid = df[~df['barcode'].str.match(r'^barcode\d{2}$')]
    if not invalid.empty:
        bad_vals = invalid['barcode'].unique().tolist()
        msg = (f"ERROR: Invalid barcode names: {bad_vals}.\n"
               f"Please ensure that all barcodes are listed as barcodeNN, e.g barcode01")
        send_error_email("MABA16S error", msg, recipients)
        sys.exit(msg)

    return df

### Function to rename the directories containing the split ONT fastq files
def renamer(input_dir, df):
    for _, row in df.iterrows():
        sample_name = str(row['sample']).strip()
        barcode_name = str(row['barcode']).strip()

        src = os.path.join(input_dir, barcode_name)
        dst = os.path.join(input_dir, sample_name)

        if not os.path.exists(src):
            print(f"WARNING: Folder '{src}' not found. Skipping.")
            continue

        if os.path.exists(dst):
            print(f"WARNING: Destination '{dst}' already exists. Skipping.")
            continue

        shutil.move(src, dst)
        print(f"Renamed {barcode_name} → {sample_name}")


def main():
    
    # Define arguments
    parser = argparse.ArgumentParser(description="Rename ONT barcode folders using a sample sheet.")

    parser.add_argument("-i", "--input-dir", required=True,
                        help="Directory containing barcodeNN folders")

    parser.add_argument("-s", "--samplesheet", required=True,
                        help="Excel file with columns 'sample' and 'barcode'")

    parser.add_argument("--technician",
                        help="Email address of the technician to notify on error")

    parser.add_argument("--mmb", 
                        help="Email address of the MMB contact to notify on error")

    args = parser.parse_args()

    # Carry out functionality
    df = validate_samplesheet(samplesheet=args.samplesheet, technician=args.technician, mmb=args.mmb)

    renamer(input_dir=args.input_dir, df=df)

    return 0


if __name__ == "__main__":
    sys.exit(main())

