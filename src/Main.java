public class Main {
    public static void main(String[] args) {
        ProbabilisticSequence seq = new ProbabilisticSequence("src/files/BoreoEutherian_chr22.txt",
                "src/files/BoreoEutherian_chr22_prob.txt", ScoringScheme.POWER1_5);
        //ProbabilisticSequence seq = new ProbabilisticSequence("src/files/test1_seq.txt",
        //        "src/files/test1_prob.txt");
        ProbabilisticBLASTTester tester = new ProbabilisticBLASTTester(seq, 200, 11, 1000, 0.07);
        tester.run();
    }
}