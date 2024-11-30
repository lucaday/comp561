public class ProbabilisticAlignment {

    private final String S;

    private final ProbabilisticSequence T;

    private final double score;

    private int genomeIndex;

    private double E;

    public ProbabilisticAlignment(String S, ProbabilisticSequence T, double score) {
        this.S = S;
        this.T = T;
        this.score = score;
    }

    public void setE(double K, int genomeLength, int qLength, double L) {
        E = K * genomeLength * qLength * Math.exp(-L * score);
    }

    public void setGenomeIndex(int genomeIndex) {
        this.genomeIndex = genomeIndex;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append("Alignment with E value: ").append(E).append(" at index ").append(genomeIndex).append("\n\t\t");
        for (char c : S.toCharArray()) {
            s.append(c).append("\t\t");
        }
        s.append("\n");
        s.append(T.toString());
        return s.toString();
    }

    public int getGenomeIndex() {
        return genomeIndex;
    }

    public double getEScore() {
        return E;
    }
}