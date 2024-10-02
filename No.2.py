def amino_acid_to_codons(amino_acids):
    """Return a list of possible codons for each amino acid."""
    codon_table = {
        'A': ['GCU', 'GCC', 'GCA', 'GCG'],
        'C': ['UGU', 'UGC'],
        'D': ['GAU', 'GAC'],
        'E': ['GAA', 'GAG'],
        'F': ['UUU', 'UUC'],
        'G': ['GGU', 'GGC', 'GGA', 'GGG'],
        'H': ['CAU', 'CAC'],
        'I': ['AUU', 'AUC', 'AUA'],
        'K': ['AAA', 'AAG'],
        'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
        'M': ['AUG'],
        'N': ['AAU', 'AAC'],
        'P': ['CCU', 'CCC', 'CCA', 'CCG'],
        'Q': ['CAA', 'CAG'],
        'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'S': ['UCU', 'UCC', 'UCA', 'UCG'],
        'T': ['ACU', 'ACC', 'ACA', 'ACG'],
        'V': ['GUU', 'GUC', 'GUA', 'GUG'],
        'W': ['UGG'],
        'Y': ['UAU', 'UAC'],
    }

    codon_combinations = []
    for amino_acid in amino_acids:
        if amino_acid in codon_table:
            codon_combinations.append(codon_table[amino_acid])
        else:
            return []

    # Generate all combinations of codons for the amino acids manually
    def cartesian_product(arrays):
        if not arrays:
            return [[]]
        result = []
        for item in arrays[0]:
            for prod in cartesian_product(arrays[1:]):
                result.append([item] + prod)
        return result

    # Get the product of the codon lists
    combinations = cartesian_product(codon_combinations)

    # Join the combinations into strings
    return [''.join(combination) for combination in combinations]

def count_codons(mrna):
    """Count the frequency of each codon in the mRNA string."""
    codon_count = {}
    for i in range(0, len(mrna), 3):
        codon = mrna[i:i+3]
        if codon in codon_count:
            codon_count[codon] += 1
        else:
            codon_count[codon] = 1
    return codon_count

def main():
    # Get user input for amino acids
    input_amino_acids = input("Input Aminoacid = ").strip().upper()

    if len(input_amino_acids) > 3 or not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in input_amino_acids):
        print("Invalid input! Please enter up to 3 valid amino acids.")
        return

    # Get all possible mRNA combinations
    codon_combinations = amino_acid_to_codons(input_amino_acids)

    # Count and print codon frequencies
    for mrna in codon_combinations:
        codon_count = count_codons(mrna)
        print(f"mRNA = {mrna}")
        for codon, count in codon_count.items():
            print(f"{codon} = {count}")

# Run the main function
if __name__ == "__main__":
    main()

