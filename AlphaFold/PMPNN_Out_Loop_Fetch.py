# Raw data provided by the user
raw_data = """
DPLSNNGKPG	580	0.2	0.601
DPLSNNGLPG	179	0.2	0.606
DPLTNNGKPG	107	0.2	0.634
DPLTNNGLPG	38	0.2	0.639
DPLSNNGRPG	20	0.2	0.666
DPLSNKGKPG	13	0.2	0.679
DPLSNNGQPG	13	0.2	0.673
DPLSNNGNPG	12	0.2	0.696
DPLSDNGKPG	5	0.2	0.679
DPLTNNGQPG	3	0.2	0.692
DPLSNNGKRG	3	0.2	0.725
DPLSNNGKLG	3	0.2	0.725
DPLSNNGEPG	3	0.2	0.709
DPLTNKGKPG	3	0.2	0.717
DPLYNNGKPG	2	0.2	0.72
DPLSNKGLPG	2	0.2	0.691
DPLSNRGKPG	2	0.2	0.718
DPLTNNGRPG	2	0.2	0.717
DMLSNFGYPG	1	0.2	1.548
DPLYNRGKPG	1	0.2	0.824
DPLYNNGLPG	1	0.2	0.74
DPLTNNGSPG	1	0.2	0.817
DPLWNNGQPG	1	0.2	0.79
DPLSNRGLPG	1	0.2	0.718
DPLSNNGHPG	1	0.2	0.716
DPLTNKGLPG	1	0.2	0.714
DPLSNNGKVG	1	0.2	0.809
DPLWNNGKPG	1	0.2	0.708
DPLTNNGHPG	1	0.2	0.76

"""

data_with_scores = []

# Process the raw data string
for line in raw_data.strip().split('\n'):
    parts = line.split()
    if parts:
        sequence = parts[0]
        score = float(parts[3])
        data_with_scores.append((sequence, score))

# Sort the list by score (the second element of the tuple) in ascending order
data_with_scores.sort(key=lambda x: x[1])

# Get the top 10 sequences
top_10_sequences = [item[0] for item in data_with_scores[:10]]

print(top_10_sequences)