package constraint;

import java.util.*;

/**
 * ConstrainedKaos  18.06.20
 * This class contains either constrained or allowed codewords in cgr matrix
 *
 * @author Hannah Franziska Löchel
 */

public class Constraints {

    /*
     *  constraints as an mCGR object
     */
    private mCGR result;
    /*
     * codeword length
     */
    private int wordlength;

    /**
     * Constructor for empty constraint object
     *
     * @param length of the desired codeword
     */
    public Constraints(int length) {
        this.wordlength = length;
        this.result = new mCGR();
        this.result.setLen(this.wordlength);
    }


    /**
     * Constructor for one constraint
     *
     * @param length    of the desired codeword
     * @param constrain String with the forbidden subsequence
     */
    public Constraints(int length, String constrain) {
        mCGR mCGR = new mCGR(constrain);
        this.wordlength = length;
        createMatrix(constrain.length(), mCGR);
    }

    /**
     * Constructor for constrained homopolymeres in a sparsematrix and constraints as fasta file
     *
     * @param length of the desired codewords
     * @param input  sequences as HashMap, with length as key and sequences as values
     */

    public Constraints(int length, HashMap<Integer, HashSet<String>> input) {

        ArrayList<Integer> arrayList = new ArrayList<>(input.keySet());
        Collections.sort(arrayList);
        int zoom = arrayList.get(0);
        mCGR mCGR = new mCGR(input.get(zoom), zoom);

        this.wordlength = length;
        while (zoom < this.wordlength) {

            if (arrayList.contains(zoom)) {
                mCGR mCGR2 = new mCGR(input.get(zoom), zoom);
                mCGR = mCGR.add(mCGR2);


            }
            mCGR nextMatrix = mCGR.doubleSize();
            mCGR = nextMatrix.tiling(mCGR);


            zoom++;
        }

        if(arrayList.contains(zoom) && arrayList.size()>1){
            mCGR mCGR2 = new mCGR(input.get(zoom), zoom);
            mCGR = mCGR.add(mCGR2);
        }
        this.result = mCGR;
    }

    /**
     * Constructor for constrained homopolymeres in a sparsematrix
     *
     * @param length of the desired codewords
     * @param hp     length of the forbidden homopolymers
     */
    public Constraints(int length, int hp) {
        mCGR mCGR = new mCGR(hp);
        this.wordlength = length;
        createMatrix(hp, mCGR);
    }


    /**
     * Construction of the matrix within constraints
     *
     * @param zoom represents the length of words in the mCGR
     * @param mCGR for tiling
     */
    private void createMatrix(int zoom, mCGR mCGR) {
        while (zoom < this.wordlength) {


            mCGR nextMatrix = mCGR.doubleSize();
            mCGR = nextMatrix.tiling(mCGR);

            zoom++;
        }
        this.result = mCGR;

    }


    /**
     * Filter for the allowed sequences within GC constraint
     *
     * @param gcContent GCContend object
     */
    public void filterGC(GCContent gcContent) {
        this.result = this.getMatrix().filter(gcContent);

    }


    /**
     * Method to calculate the ratio of words in the sparsematrix to all words of this length
     *
     * @return ratio of codewords to possible combinations of this word length
     */
    public double ratio() {
        int all = this.getMatrix().getLen() * this.getMatrix().getLen();
        int allowed = 0;
        HashMap<Integer, HashSet<Integer>> map = this.getMatrix().getRows();
        for (int element : map.keySet()) {
            allowed = allowed + map.get(element).size();
        }
        System.out.println("allowed: " + allowed + " from: " + all);

        return (double) allowed / all;
    }

    /**
     * Method to get the codeword length of an constraint object
     *
     * @return length of codewords
     */
    public int getWordlength() {
        return this.wordlength;
    }


    /**
     * Method to get the object containing the results
     *
     * @return mCGR Object with result of constrained / allowed sequences
     */
    public mCGR getMatrix() {
        return this.result;
    }

    /**
     * Method to delete sequences, which can form an undesired motif
     *
     * @param links HashSet containing mCGRs with positions forming undesired motifs when concatenated
     */
    public void filterMotifs(Link links) {
        this.result = this.getMatrix().filter(links);
    }
}