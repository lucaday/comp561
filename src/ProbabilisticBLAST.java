import java.util.List;
import java.util.Map;

public class ProbabilisticBLAST implements Runnable {

    private final String q;

    private final ProbabilisticSequence genome;

    private final Map<String, List<Integer>> indices;

    private final int wordLength;

    private final double delta;

    private final double threshold;

    private final double K;

    private final double L;

    public ProbabilisticBLAST(String q, ProbabilisticSequence genome, int wordLength) {
        this.q = q;
        this.genome = genome;
        this.indices = genome.buildIndices(wordLength);
        this.wordLength = wordLength;
        this.delta = 10;
        this.threshold = 1;
        this.K = 1;
        this.L = 1;
    }

    @Override
    public void run() {
        // For each w-mer in q
        List<Integer> genomeIndices;
        for (int i = 0; i < q.length() - wordLength; i++) {
            genomeIndices = this.indices.get(q.substring(i, i + wordLength));
            if (genomeIndices == null) continue;

            // Perform ungapped extension at each genomeIndex where this substring occurs ('seed' is determined by i and index)
            for (Integer genomeIndex : genomeIndices) {
                if (ungappedExtension(i, genomeIndex)) {
                    System.out.println("Found HSP at " + genomeIndex);
                }
            }
        }
    }

    private boolean ungappedExtension(int qStartIndex, int genomeStartIndex) {
        // Ungapped extension from the right
        int qIndex = qStartIndex;
        int genomeIndex = genomeStartIndex;
        double currScore = 0.0;
        double maxScoreRight = 0;
        while (qIndex < q.length() && genomeIndex < genome.length()) {
            currScore += scoreChar(q.charAt(qIndex), genomeIndex);
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
            currScore += scoreChar(q.charAt(qIndex), genomeIndex);
            if (currScore >= maxScoreLeft) maxScoreLeft = currScore; // Greater than OR EQUAL TO! Because a probability of 0.5 gives a score of 0
            if (currScore < maxScoreLeft - delta) break;
            qIndex--;
            genomeIndex--;
        }

        // Get total score from max left and right scores, return true if E is less than threshold (HSP found!)
        double totalScore = maxScoreLeft + maxScoreRight;
        // E = K m n e ^ (-L * s), s = alignment score, m = |D|, n = |Q|, L and K depend on scoring scheme
        double E = K * genome.length() * q.length() * Math.exp(-L * totalScore);

        return E < threshold;
    }

    private double scoreChar(char c, int genomeIndex) {
        // using 2x - 1 to linearly rescale probabilities from [0, 1] to [-1, 1]
        return 2 * genome.getProbabilities().get(c).get(genomeIndex) - 1;
    }
}
