import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.stream.Stream;

public class ProbabilisticSequence {

    private final Map<Character, List<Double>> probabilities;

    private final ScoringScheme scoringScheme;

    private ProbabilisticSequence(Map<Character, List<Double>> probabilities, ScoringScheme scoringScheme) {
        this.probabilities = probabilities;
        this.scoringScheme = scoringScheme;
    }

    public ProbabilisticSequence() {
        this.probabilities = new HashMap<>();
        for (char nucleotide : List.of('A', 'C', 'G', 'T')) {
            this.probabilities.put(nucleotide, new ArrayList<>());
        }
        this.scoringScheme = null;
    }

    /**
     * Constructs a `ProbabilisticSequence` object by loading sequence file and probability file which builds a
     * probability matrix representing each nucleotide (A, C, G, T) at every position in the sequence.
     *
     * @param sequenceFileName the name of the file containing the DNA sequence (as a single string).
     * @param probabilitiesFileName the name of the file containing space-delimited probabilities for
     *                              the most likely nucleotide at each position in the sequence (the rest of the
     *                              nucleotides will be assumed to all have equal probabilities)
     */
    public ProbabilisticSequence(String sequenceFileName, String probabilitiesFileName, ScoringScheme scoringScheme) {
        probabilities = new HashMap<>();
        this.scoringScheme = scoringScheme;
        List.of('A', 'C', 'G', 'T').forEach(n -> probabilities.put(n, new ArrayList<>()));

        try (BufferedReader sequenceReader = new BufferedReader(new FileReader(sequenceFileName));
             BufferedReader probabilityReader = new BufferedReader(new FileReader(probabilitiesFileName))) {
            String sequence = sequenceReader.readLine();
            String probabilitiesString = probabilityReader.readLine();

            if (sequence == null) {
                throw new IllegalArgumentException("Sequence file is empty");
            }
            if (probabilitiesString == null) {
                throw new IllegalArgumentException("Probabilities file is empty");
            }

            List<Double> probabilityList = Stream.of(probabilitiesString.split(" "))
                    .map(Double::parseDouble)
                    .toList();

            if (sequence.length() != probabilityList.size()) {
                throw new IllegalArgumentException("File sizes do not match");
            }

            for (int i = 0; i < sequence.length(); i++) {
                char c = sequence.charAt(i);
                double mainProb = probabilityList.get(i);
                double otherProb = (1.0 - mainProb) / 3.0;

                probabilities.get('A').add(c == 'A' ? mainProb : otherProb);
                probabilities.get('C').add(c == 'C' ? mainProb : otherProb);
                probabilities.get('G').add(c == 'G' ? mainProb : otherProb);
                probabilities.get('T').add(c == 'T' ? mainProb : otherProb);
            }

        } catch (IOException e) {
            throw new RuntimeException("Error reading files: " + e.getMessage(), e);
        }
    }

    /**
     * Builds an index mapping substrings (of a specified length) to their positions in a probabilistic sequence.
     * <p>
     * This method selects nucleotides based on which has the highest probability at a given position. It is assumed
     * there are no positions where two nucleotides have the same probability (i.e. no tiebreaks).
     *
     * @param wordLength the length of the words to extract and index.
     * @return a map where keys are substrings of length `wordLength` and values are lists of starting
     *         positions where these substrings occur in the sequence.
     */
    public Map<String, List<Integer>> buildIndices(int wordLength) {
        HashMap<String, List<Integer>> indexMap = new HashMap<>();
        StringBuilder sequence = new StringBuilder();
        for (int i = 0; i < length(); i++) {
            sequence.append(getMaxProbChar(i));
            if (sequence.length() > wordLength) sequence.deleteCharAt(0);
            if (sequence.length() == wordLength) indexMap.computeIfAbsent(sequence.toString(), k -> new ArrayList<>()).add(i - wordLength + 1);
        }
        return indexMap;
    }

    /**
     * Generates a random sequence based on the probabilistic sequence matrix, incorporating optional
     * random insertions and deletions.
     *
     * @param startingIndex  the starting index in the probabilistic sequence.
     * @param length       the length to go through the probabilistic sequence
     * @param insertionProb  the probability of inserting a random nucleotide at any position.
     * @param deletionProb   the probability of deleting a nucleotide from the sequence.
     * @param r           java util random class (for reproducibility).
     * @return a randomly generated sequence as a String.
     */
    public String generateRandomSequence(int startingIndex, int length, double insertionProb, double deletionProb, Random r) {
        StringBuilder sequence = new StringBuilder();
        for (int i = startingIndex; i < startingIndex + length; i++) {
            if (r.nextDouble() < deletionProb) continue;
            if (r.nextDouble() < insertionProb) sequence.append(List.of('A', 'C', 'G', 'T').get(r.nextInt(4)));

            double randomValue = r.nextDouble();
            double probA = probabilities.get('A').get(i);
            double probC = probabilities.get('C').get(i);
            double probG = probabilities.get('G').get(i);

            if (randomValue < probA) {
                sequence.append('A');
            } else if (randomValue < probA + probC) {
                sequence.append('C');
            } else if (randomValue < probA + probC + probG) {
                sequence.append('G');
            } else {
                sequence.append('T');
            }
        }
        return sequence.toString();
    }

    public ProbabilisticSequence subSequence(int startIndex, int endIndex) {
        if (startIndex < 0) startIndex = 0;
        if (endIndex > length()) endIndex = length();
        Map<Character, List<Double>> subProbabilities = new HashMap<>();
        for (char nucleotide : List.of('A', 'C', 'G', 'T')) {
            subProbabilities.put(nucleotide, new ArrayList<>(probabilities.get(nucleotide).subList(startIndex, endIndex)));
        }
        return new ProbabilisticSequence(subProbabilities, scoringScheme);
    }

    public double scoreChar(char c, int i) {
        switch (scoringScheme) {
            case POWER1 -> {
                return 2 * probabilities.get(c).get(i) - 1;
            }
            case POWER0_5 -> {
                return 2 * Math.pow(probabilities.get(c).get(i), 0.5) - 1;
            }
            case POWER1_5 -> {
                return 2 * Math.pow(probabilities.get(c).get(i), 1.5) - 1;
            }
            default -> throw new IllegalArgumentException("Unrecognized scoring scheme");
        }
    }

    public int length() {
        return probabilities.get('A').size();
    }

    private char getMaxProbChar(int index) {
        return Stream.of('A', 'C', 'G', 'T')
                .max((c1, c2) -> Double.compare(probabilities.get(c1).get(index), probabilities.get(c2).get(index)))
                .orElseThrow();
    }

    public void insertEmpty() {
        for (char nucleotide : List.of('A', 'C', 'G', 'T')) {
            probabilities.get(nucleotide).add(0, 0.);
        }
    }

    public void insertCopy(ProbabilisticSequence seq, int index) {
        for (char nucleotide : List.of('A', 'C', 'G', 'T')) {
            probabilities.get(nucleotide).add(0, seq.probabilities.get(nucleotide).get(index));
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (char nucleotide : List.of('A', 'C', 'G', 'T')) {
            sb.append(nucleotide).append(":\t\t");
            for (double probability : probabilities.get(nucleotide)) {
                sb.append(String.format("%.3f", probability)).append("\t");
            }
            sb.append("\n");
        }

        return sb.toString();
    }
}