from Bio import SeqIO
import itertools
import random
lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]

example = ["ACGT", "ATGT", "CCGA"]

#TASK1
def count_matrix(motifs):
    l = len(motifs[0])
    counts = {"A": [0] * l, "C": [0] * l, "G": [0] * l, "T": [0] * l}
    for motif in motifs:
        for i in range(l):
            position = motif[i]
            counts[position][i] += 1
    return counts
#print(count_matrix(example))

def score(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    total_score = 0
    for i in range(l):
        column = [counts["A"][i], counts["C"][i], counts["G"][i], counts["T"][i]]
        total_score += max(column)
    return total_score
#print(score(example))

def consensus(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])

    result_consensus = ""

    for i in range(l):
        b_base = "A"
        b_count = counts["A"][i]

        for base in "CGT":
            if counts[base][i] > b_count:
                b_base = base
                b_count = counts[base][i]

        result_consensus += b_base

    return result_consensus
#print(consensus(example))

def hamming_distance(a, b):
    distance = 0
    for i in range(len(a)):
        if a[i] != b[i]:
            distance += 1
    return distance

def total_distance(pattern, sequences):
    l = len(pattern)
    total_d = 0
    for sequence in sequences:
        best_distance = l
        for i in range(len(sequence) - l + 1):
            window = sequence[i:i + l]
            distance = hamming_distance(pattern, window)
            if distance < best_distance:
                best_distance = distance
        total_d += best_distance

    return total_d
#print(total_distance("AC", ["GACT", "TTAG"]))

#red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
#best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

#print(score(red))                                # 26
#print(consensus(best), score(best))              # AGATAG 34
#print(hamming_distance("TAAGTT", "TGAATT"))      # 2
#print(total_distance("TGCGTT", lecture_dna))     # 13

#TASK2

class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.l = len(motifs[0])
        length = len(motifs)
        counts = count_matrix(motifs)
        self.ppm = {"A" :[0.0] * self.l, "C" :[0.0] * self.l, "G" :[0.0] * self.l, "T" :[0.0] * self.l}

        denominator = length + 4 * pseudocount
        for base in "ACGT":
            for i in range(self.l):
                self.ppm[base][i] = (counts[base][i] + pseudocount) / denominator

    def lmer_probability(self,lmer):
        probability = 1.0
        for i in range(self.l):
            base = lmer[i]
            probability *= self.ppm[base][i]

        return probability

    def most_probable_lmer(self, sequence):
        best_lmer = sequence[0:self.l]
        best_probability = self.lmer_probability(best_lmer)

        for i in range(1, len(sequence) - self.l +1):
            current_lmer = sequence[i:i + self.l]
            current_probability = self.lmer_probability(current_lmer)
            if current_probability > best_probability:
                best_lmer = current_lmer
                best_probability = current_probability
        return best_lmer

    def consensus(self):
        result = ""
        for i in range(self.l):
            b_base = "A"
            b_probability = self.ppm["A"][i]

            for base in "CGT":
                if self.ppm[base][i] > b_probability:
                    b_base = base
                    b_probability = self.ppm[base][i]
            result += b_base
        return result

profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
#print(round(profile.lmer_probability("ATGCGTA"), 4))
#print(profile.consensus())

#two = MotifProfile(["GTAC", "TTAA"])
#print(two.most_probable_lmer("ACTGGATGACCC"))    # TGAC
#print(round(two.lmer_probability("TGAC"), 4))         # 0.0093

#TASK3
class MotifFinder:
    def __init__(self, sequences, l, seed=None):
        self.sequences = sequences
        self.l = l
        self.rng = random.Random(seed)
        self.windows = []
        for sequence in self.sequences:
            sequence_windows = []
            for i in range(len(sequences) - self.l +1):
                window = sequence[i:i + self.l]
                sequence_windows.append(window)
            self.windows.append(sequence_windows)

    def total_distance(self, pattern):
        total = 0
        for sequence_windows in self.windows:
            best_distance = self.l
            for window in sequence_windows:
                distance = hamming_distance(pattern,window)
                if distance < best_distance:
                    best_distance = distance
            total += best_distance
            return total

    def median_string(self):
        best_pattern = None
        best_distance = None
        for combination in itertools.product("ACGT",repeat=self.l):
            pattern = "".join(combination)
            distance = self.total_distance(pattern)

            if best_distance is None or distance < best_distance:
                best_pattern = pattern
                best_distance = distance

        return best_pattern, best_distance

    def randomized_search(self):
        current_motifs = []
        for sequnce_windows in self.windows:
            random_lmer = self.rng.choice(sequnce_windows)
            current_motifs.append(random_lmer)
        current_score = score(current_motifs)
        while True:
            profile = MotifProfile(current_motifs,pseudocount=1)
            new_motifs = []

            for sequence in self.sequences:
                best_lmer = profile.most_probable_lmer(sequence)
                new_motifs.append(best_lmer)
            new_score = score(new_motifs)

            if new_score > current_score:

                current_motifs = new_motifs
                current_score = new_score

            else:
                return current_motifs, current_score

    def best_of(self, runs):

        best_motifs = None
        best_score = -1
        for _ in range(runs):
            current_motifs, current_score = (self.randomized_search())

            if current_score > best_score:
                best_motifs = current_motifs
                best_score = current_score

        return best_motifs, best_score




