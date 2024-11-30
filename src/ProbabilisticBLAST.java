import java.util.*;

public class ProbabilisticBLAST {

    private final String q;

    private final ProbabilisticSequence genome;

    private final Map<String, List<Integer>> indices;

    private final int wordLength;

    private final double delta;

    private final double threshold;

    private final double K;

    private final double L;

    public ProbabilisticBLAST(String q,
                              ProbabilisticSequence genome,
                              Map<String, List<Integer>> indices,
                              int wordLength,
                              double threshold,
                              double K,
                              double L) {
        this.q = q;
        this.genome = genome;
        this.indices = indices;
        this.wordLength = wordLength;
        this.delta = 10;
        this.threshold = threshold;
        this.K = K;
        this.L = L;
    }

    public List<ProbabilisticAlignment> getAlignments() {
        // For each w-mer in q
        List<Integer> genomeIndices;
        Set<Integer> exploredIndices = new HashSet<>();
        List<ProbabilisticAlignment> alignments = new ArrayList<>();
        for (int i = 0; i < q.length() - wordLength; i++) {
            genomeIndices = this.indices.get(q.substring(i, i + wordLength));
            if (genomeIndices == null) continue;

            // Perform ungapped extension at each genomeIndex where this substring occurs ('seed' is determined by i and index)
            for (Integer genomeIndex : genomeIndices) {
                double E = ungappedExtension(i, genomeIndex);
                if (E < threshold) {
                    int indexToExplore = genomeIndex - i;
                    if (exploredIndices.add(indexToExplore)) {
                        ProbabilisticNeedlemanWunsch needlemanWunsch = new ProbabilisticNeedlemanWunsch(q, genome.subSequence(indexToExplore, indexToExplore + q.length()));
                        ProbabilisticAlignment alignment = needlemanWunsch.getAlignment();
                        alignment.setE(K, genome.length(), q.length(), L);
                        alignment.setGenomeIndex(indexToExplore);
                        if (alignment.getEScore() < threshold) {
                            alignments.add(alignment);
                        }
                    }
                }
            }
        }
        return alignments;
    }

    private double ungappedExtension(int qStartIndex, int genomeStartIndex) {
        // Ungapped extension from the right
        int qIndex = qStartIndex;
        int genomeIndex = genomeStartIndex;
        double currScore = 0.0;
        double maxScoreRight = 0;
        while (qIndex < q.length() && genomeIndex < genome.length()) {
            currScore += genome.scoreChar(q.charAt(qIndex), genomeIndex);
            if (currScore >= maxScoreRight) maxScoreRight = currScore; // Greater than OR EQUAL TO! Because a probability of 0.5 gives a score of 0
            if (currScore < maxScoreRight - delta) break;
            qIndex++;
            genomeIndex++;
        }

        // Ungapped extension from the left
        qIndex = qStartIndex - 1;
        genomeIndex = genomeStartIndex - 1;
        currScore = 0.0;
        double maxScoreLeft = 0;
        while (qIndex >= 0 && genomeIndex >= 0) {
            currScore += genome.scoreChar(q.charAt(qIndex), genomeIndex);
            if (currScore >= maxScoreLeft) maxScoreLeft = currScore; // Greater than OR EQUAL TO! Because a probability of 0.5 gives a score of 0
            if (currScore < maxScoreLeft - delta) break;
            qIndex--;
            genomeIndex--;
        }

        // Get total score from max left and right scores, return true if E is less than threshold (HSP found!)
        double totalScore = maxScoreLeft + maxScoreRight;
        // E = K m n e ^ (-L * s), s = alignment score, m = |D|, n = |Q|, L and K depend on scoring scheme
        return K * genome.length() * q.length() * Math.exp(-L * totalScore);
    }
}
