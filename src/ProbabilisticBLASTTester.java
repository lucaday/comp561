import java.util.List;
import java.util.Map;
import java.util.Random;

public class ProbabilisticBLASTTester implements Runnable {

    private final ProbabilisticSequence genome;

    private final int qLength;

    private final int wordLength;

    private final int nTests;

    private final double gapProb;

    public ProbabilisticBLASTTester(ProbabilisticSequence genome,
                                    int qLength,
                                    int wordLength,
                                    int nTests,
                                    double gapProb) {
        this.genome = genome;
        this.qLength = qLength;
        this.wordLength = wordLength;
        this.nTests = nTests;
        this.gapProb = gapProb;
    }

    @Override
    public void run() {
        int max = genome.length() - qLength;
        Random random = new Random(1);
        String q;
        ProbabilisticBLAST blast;
        Map<String, List<Integer>> indices = genome.buildIndices(wordLength);
        int startIndex;
        double pos = 0.;
        double falsePos = 0.;
        for (int i = 0; i < nTests; i++) {
            startIndex = random.nextInt(0, max);
            q = genome.generateRandomSequence(startIndex, qLength, gapProb, gapProb, random);
            blast = new ProbabilisticBLAST(q, genome, indices, wordLength, 1, 1, 1);
            List<ProbabilisticAlignment> alignments = blast.getAlignments();
            boolean posTrue = false;
            for(ProbabilisticAlignment alignment : alignments) {
                if (alignment.getGenomeIndex() >= (startIndex - qLength * gapProb) &&
                        alignment.getGenomeIndex() <= (startIndex + qLength * gapProb)) {
                    if (!posTrue) {
                        posTrue = true;
                        pos++;
                    }
                } else {
                    falsePos++;
                }
            }
        }

        System.out.println(pos / nTests);
        System.out.println(falsePos / nTests);
    }
}
