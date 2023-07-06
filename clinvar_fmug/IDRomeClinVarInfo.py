# Copyright (C) 2023 Nicolas Jonsson <nicolas.jonsson@bio.ku.dk>

"""Module for data handling in the IDRome project

This module implements classes for importing local application 
modules (LocalAppImports), parsing and processing data files 
(DataProcessor) in the context of the IDRome project. 

The DataProcessor class handles loading data from a CSV file, processing 
ClinVar benign and pathogenic annotations in the IDR regions, and exporting 
the processed data to a CSV file. 

The LocalAppImports class handles logging configuration and imports specific 
PRISM modules from a specified directory. The imported modules are then assigned 
as instance methods.

In general, this script is used to process IDR fragment data with relevant ClinVar 
information and output a structured dataset for further analysis in the IDRome project.

"""

import glob
import logging as log
import os
import sys
import warnings
import pandas as pd
import tqdm
# Ignore all warnings
warnings.filterwarnings("ignore")


class LocalAppImports:
    """A class to handle local application imports and logging configuration"""

    def __init__(self, log_message: str):
        self.log_message = log_message
        self.configure_logging()
        self.import_prism_modules()

    def configure_logging(self):
        """Configure logging based on the provided log_message"""
        log_levels = {
            "verbose": log.INFO,
            "debug": log.WARNING,
            "error": log.ERROR,
        }

        log.basicConfig(
            format="%(levelname)s:%(message)s",
            datefmt="%m/%d/%Y %I:%M:%S %p",
            level=log_levels.get(self.log_message, log.ERROR),
        )

    def import_prism_modules(self):
        """Import PRISM modules and assign them as instance methods"""
        logger = log.getLogger(__name__)
        basepath = "/storage1/jonsson/software/"

        try:
            sys.path.insert(1, os.path.join(basepath, "prism/scripts/"))
            from PrismData import PrismParser, VariantData
            sys.path.insert(1, os.path.join(basepath, "KULL_PRISM/software/scripts"))
            from prism_parser_helper import (
                read_from_prism,
                read_prism,
                merge_prism,
                merge_prism_df,
                write_prism,
            )
        except (ModuleNotFoundError, ImportError) as e:
            logger.error(f"{type(e)} failure")
            print(e)
        else:
            logger.info("Import succeeded")

        self.read_from_prism = read_from_prism
        self.read_prism = read_prism
        self.merge_prism = merge_prism
        self.merge_prism_df = merge_prism_df
        self.write_prism = write_prism


class DataProcessor:
    """A class to handle loading, processing, and exporting data"""

    def __init__(self, local_app_imports: LocalAppImports):
        self.local_app_imports = local_app_imports
        self.dataframe = None

    def load_dataframe(self, filepath: str):
        """Load dataframe from a csv file and rename the first column"""
        self.dataframe = pd.read_csv(filepath)
        self.dataframe.rename(columns={"Unnamed: 0": "IDRFragments"}, inplace=True)

    def process_data(self) -> list:
        """
        Process the loaded data to extract ClinVar benign and pathogenic 
        annotations in the IDR regions and return a list of processed data
        """
        prismDir = "/storage1/shared/data/prism"
        db_columns = [
            "IDRFragments",
            "UniProt_ID",
            "N",
            "fasta",
            "N_term",
            "C_term",
            "first",
            "last",
            "protein_name",
            "gene_name",
        ]
        db_values = {
            column: self.dataframe[column].to_numpy() for column in db_columns
        }
        ClinVarInfo = []
        
        # Loop through IDRFragments with a progress bar
        for i, frag in tqdm.tqdm(
            enumerate(db_values["IDRFragments"]), total=len(db_values["IDRFragments"])
        ):
            UniProt_ID, first, last = frag.split("_")
            first, last = int(first), int(last)
            fasta = db_values["fasta"][i]
            
            ClinVarFiles = glob.glob(
                os.path.join(
                    prismDir,
                    UniProt_ID[0:2],
                    UniProt_ID[2:4],
                    UniProt_ID[4:6],
                    f"prism_clinvar_*_{UniProt_ID}*.txt",
                )
            )

            # Process ClinVar files if they exist
            for file in ClinVarFiles:
                meta_data, dataframe = self.local_app_imports.read_from_prism(file)
                dataframe = dataframe[dataframe["Clinvar_signifiance"].isin(["benign", "pathogenic"])]

                dataframe["aa_ref"] = dataframe["aa_ref"].explode()
                dataframe["resi"] = dataframe["resi"].explode()
                dataframe["aa_var"] = dataframe["aa_var"].explode()

                if fasta in meta_data["protein"]["sequence"]:
                    # Filtering and sorting data
                    dataframe = dataframe[
                        (dataframe["resi"] >= first)
                        & (dataframe["resi"] <= last)
                        & (dataframe["aa_var"] != "*")
                    ]

                    # Data mining for ClinVarInfo
                    clinData = self._prepare_clinData(frag, UniProt_ID, first, last, fasta, dataframe)
                    ClinVarInfo.append(clinData)

            # If no ClinVar files exist for the fragment
            if not ClinVarFiles:
                emptyData = [[], [], [], []]
                clinData = [frag, UniProt_ID, first, last, fasta] + emptyData
                ClinVarInfo.append(clinData)

        return ClinVarInfo

    def _prepare_clinData(self, frag: str, UniProt_ID: str, first: int, last: int, fasta: str, dataframe: pd.DataFrame) -> list:
        """Prepare a list of ClinVar data for a given fragment"""
        if dataframe.empty:
            emptyData = [[], [], [], []]
            clinData = [frag, UniProt_ID, first, last, fasta] + emptyData
        else:
            dataframe.sort_values(by=["resi"], inplace=True)

            # Extracting values from filtered temporary dataframe
            variant = dataframe["variant"].tolist()
            ClinVarID = dataframe["Clinvar_ID"].tolist()
            ClinVarSignificance = dataframe["Clinvar_signifiance"].tolist()
            ClinVarReviewStatus = dataframe["Clinvar_review_status"].tolist()

            clinData = [
                frag,
                UniProt_ID,
                first,
                last,
                fasta,
                variant,
                ClinVarID,
                ClinVarSignificance,
                ClinVarReviewStatus,
            ]
        return clinData

    def export_data(self, data: list, filename: str):
        """Export the processed data to a csv file"""
        dataframe = pd.DataFrame(
            data,
            columns=[
                "IDR_fragment",
                "UniProt_ID",
                "first",
                "last",
                "fasta",
                "variants",
                "Clinvar_ID",
                "Clinvar_signifiance",
                "Clinvar_review_status",
            ],
        )
        dataframe.to_csv(filename, index=False)


def main():
    local_app_imports = LocalAppImports(log_message="verbose")
    data_processor = DataProcessor(local_app_imports)
    data_processor.load_dataframe( ) #Insert path to the 'IDRome_DB.csv' file
    processed_data = data_processor.process_data()
    data_processor.export_data(data=processed_data, filename="./IDRome_DB_ClinVarInfo.csv")

if __name__ == "__main__":
    main()