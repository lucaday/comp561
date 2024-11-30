public class ProbabilisticNeedlemanWunsch {

    private final String A;

    private final ProbabilisticSequence B;

    public ProbabilisticNeedlemanWunsch(String A, ProbabilisticSequence B) {
        this.A = A;
        this.B = B;
    }

    public ProbabilisticAlignment getAlignment() {
        char[] a = ('-' + A).toCharArray();
        B.insertEmpty();
        Element[][] H = scoringMatrix(a, B, 2 * 0.07 - 1);
        return traceback(a, B, H);
    }

    private record Element(double value, int direction) {}

    private Element[][] scoringMatrix(char[] A, ProbabilisticSequence B, double gapCost) {
        Element[][] H = new Element[A.length][B.length()];

        H[0][0] = new Element(0, 0);
        for (int i = 1; i < A.length; i++) {
            H[i][0] = new Element(H[i-1][0].value + gapCost, 3);
        }
        for (int j = 1; j < B.length(); j++) {
            H[0][j] = new Element(H[0][j-1].value + gapCost, 1);
        }

        for (int i = 1; i < A.length; i++) {
            for (int j = 1; j < B.length(); j++) {
                double diagonal = H[i-1][j-1].value + B.scoreChar(A[i], j);
                double top = H[i-1][j].value + gapCost;
                double left = H[i][j-1].value + gapCost; // left

                double max = diagonal;
                int direction = 2; // Assume diagonal is max initially

                if (top > max) {
                    max = top;
                    direction = 3; // Top
                }

                if (left > max) {
                    max = left;
                    direction = 1; // Left
                }

                H[i][j] = new Element(max, direction);
            }
        }

        return H;
    }

    private ProbabilisticAlignment traceback(char[] A, ProbabilisticSequence B, Element[][] H) {
        StringBuilder S = new StringBuilder();
        ProbabilisticSequence T = new ProbabilisticSequence();

        int i = A.length - 1;
        int j = B.length() - 1;
        Element element = H[i][j];

        while (element.direction != 0) {
            switch (element.direction) {
                case 1 -> { // left
                    S.insert(0, '-');
                    T.insertCopy(B, j);
                    j--;
                }
                case 2 -> { // diagonal
                    S.insert(0, A[i]);
                    T.insertCopy(B, j);
                    i--;
                    j--;
                }
                case 3 -> { // up
                    S.insert(0, A[i]);
                    T.insertEmpty();
                    i--;
                }
                default -> throw new RuntimeException("Invalid direction");
            }
            element = H[i][j];
        }

        return new ProbabilisticAlignment(S.toString(), T, H[A.length - 1][B.length() - 1].value);
    }
}
