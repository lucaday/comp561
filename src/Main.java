public class Main {
    public static void main(String[] args) {
        ProbabilisticSequence seq = new ProbabilisticSequence("src/files/BoreoEutherian_chr22.txt",
                "src/files/BoreoEutherian_chr22_prob.txt");
        String q = seq.generateRandomSequence(234503, 50, 0.07, 0.07, 1);
        System.out.println(q);
        ProbabilisticBLAST blast = new ProbabilisticBLAST(q, seq, 11);
        blast.run();
    }
}