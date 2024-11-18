
import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from io import StringIO

# Streamlit App Title
st.title("CRISPR Spacer Design Tool2")

# Upload the GenBank file
genbank_file = st.file_uploader("Upload a GenBank file (.gb or .gbff)", type=["gb", "gbff"])

# User inputs for parameters
target_genes = st.text_input("Target Genes (comma-separated):", value="edd")
pam_sequence = st.text_input("PAM Sequence:", value="CC")
spacer_length = st.slider("Spacer Length:", min_value=16, max_value=50, value=32)
left_flank = st.text_input("Left Flank Sequence:", value="aggtcTcaaaac")
right_flank = st.text_input("Right Flank Sequence:", value="gtttttGAGACCa")

# Function to extract gene sequences from GenBank
def get_gene_sequences(gb_content, gene_names):
    gene_sequences = {}
    gb_file = StringIO(gb_content)
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene" and "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                if gene_name in gene_names:
                    gene_seq = feature.location.extract(record).seq
                    gene_sequences[gene_name] = gene_seq
    return gene_sequences

# Function to generate spacers with PAM
def generate_spacers_with_pam(gene_sequences, pam_seq, spacer_length):
    spacers = {}
    for gene, sequence in gene_sequences.items():
        gene_spacers = []
        for i in range(len(sequence) - spacer_length):
            if sequence[i:i + len(pam_seq)] == pam_seq:
                spacer = sequence[i + len(pam_seq):i + len(pam_seq) + spacer_length]
                gene_spacers.append(str(spacer))
            if len(gene_spacers) == 3:  # Limit to 3 spacers per gene
                break
        spacers[gene] = gene_spacers
    return spacers

# Run analysis when the button is clicked
if st.button("Generate Spacers"):
    if genbank_file is None:
        st.error("Please upload a GenBank file.")
    else:
        try:
            # Read and parse the GenBank file
            gb_content = genbank_file.getvalue().decode("utf-8")
            target_genes_list = [gene.strip() for gene in target_genes.split(",")]

            # Extract gene sequences
            gene_sequences = get_gene_sequences(gb_content, target_genes_list)

            # Check if any gene sequences were found
            if not gene_sequences:
                st.error("No target genes found in the GenBank file. Please check the gene names.")
            else:
                # Generate spacers
                spacers = generate_spacers_with_pam(gene_sequences, pam_sequence, spacer_length)

                # Prepare data for display and download
                spacer_data = []
                for gene, spacer_list in spacers.items():
                    for i, spacer in enumerate(spacer_list):
                        # Create sense oligo
                        full_spacer = f"{left_flank}{spacer}{right_flank}"
                        spacer_name_sense = f"{gene}_sp.{i + 1}_Sense"

                        # Create antisense oligo (reverse complement)
                        reverse_complement = str(Seq(full_spacer).reverse_complement())
                        spacer_name_antisense = f"{gene}_sp.{i + 1}_AntiSense"

                        # Append to spacer data
                        spacer_data.append([spacer_name_sense, full_spacer])
                        spacer_data.append([spacer_name_antisense, reverse_complement])

                # Display results
                df = pd.DataFrame(spacer_data, columns=["Name", "Sequence"])
                st.write("Generated Spacers:")
                st.dataframe(df)

                # Provide download option
                st.download_button(
                    label="Download Results as Excel",
                    data=df.to_excel(index=False, engine="openpyxl"),
                    file_name="spacers_output.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                )

        except Exception as e:
            st.error(f"An error occurred: {e}")
