import pandas as pd
import sys
from datetime import datetime
import os
from pathlib import Path

def process_single_excel_file(excel_file_path, output_dir):
    """
    Reads an Excel file, extracts the first 5 rows of each sheet,
    and saves them to a single time-stamped text file.
    Each sheet's content is separated by its name and a line separator.
    """
    try:
        # Generate time-stamped filename based on original Excel file name
        base_excel_name = os.path.splitext(os.path.basename(excel_file_path))[0]
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_filename = os.path.join(output_dir, f"{base_excel_name}_report_{timestamp}.txt")

        with open(output_filename, "w", encoding="utf-8") as f_out:
            f_out.write(f"--- Report for Excel file: {excel_file_path} ---\n")
            f_out.write(f"Report Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f_out.write(f"\n")

            excel_file = pd.ExcelFile(excel_file_path)

            for sheet_name in excel_file.sheet_names:
                f_out.write(f"{'='*60}\n")
                f_out.write(f"Sheet Name: {sheet_name}\n")
                f_out.write(f"{'='*60}\n")
                f_out.write(f"\n")

                # Read the first 5 rows of the current sheet
                df = pd.read_excel(excel_file, sheet_name=sheet_name, nrows=5)

                f_out.write("--- Column Headers ---\n")
                if not df.columns.empty:
                    for col in df.columns:
                        f_out.write(f"- {col}\n")
                else:
                    f_out.write("No columns found (sheet might be empty or malformed).\n")

                f_out.write("\n--- First 5 Rows (including header) ---\n")
                # Use to_string for better formatting in a text file
                f_out.write(df.to_string(index=False) + "\n")

                f_out.write("\n\n") # Add extra newlines for separation between sheets

        print(f"Processed '{os.path.basename(excel_file_path)}'. Report saved to: {output_filename}")

    except FileNotFoundError:
        print(f"Error: File not found at {excel_file_path}")
    except Exception as e:
        print(f"An error occurred while processing '{excel_file_path}': {e}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python header_check.py <path_to_directory_with_excel_files>")
        sys.exit(1)

    input_directory = sys.argv[1]
    # Output reports next to this script's directory
    base_dir = Path(__file__).resolve().parent
    output_base_dir = str(base_dir)

    # Create output directory if it doesn't exist
    if not os.path.exists(output_base_dir):
        os.makedirs(output_base_dir)

    print(f"--- Processing Excel files in directory: {input_directory} ---")

    if not os.path.isdir(input_directory):
        print(f"Error: '{input_directory}' is not a valid directory.")
        sys.exit(1)

    processed_count = 0
    for filename in os.listdir(input_directory):
        if filename.endswith(".xlsx"):
            full_file_path = os.path.join(input_directory, filename)
            if os.path.isfile(full_file_path):
                process_single_excel_file(full_file_path, output_base_dir)
                processed_count += 1

    print(f"\n--- Finished processing. {processed_count} Excel file(s) processed. ---")
    print(f"Reports saved in the '{output_base_dir}' directory.")