package cgr;


import constraint.Constraints;
import constraint.mCGR;
import org.apache.commons.math3.fraction.BigFraction;

import java.io.*;
import java.math.BigInteger;
import java.util.HashMap;

/**
 * ConstrainedKaos 05.11.19
 * Calculating Sequences from last CGR coordinates to DNA string
 *
 * @author Hannah Franziska Löchel
 */
public class ReverseKaos extends CGR {
    private HashMap<String, int[]> sequences;


    /**
     * Constructor for codewords from mCGR with constrains
     * @param constraints which will be translated to DNA codewords
     */
    public ReverseKaos(Constraints constraints) {
        super();
        this.sequences = new HashMap<>();
        System.out.println("Starting to translate to DNA.");
        mCGR matrix = constraints.getMatrix();
        int len = matrix.getLen();
        for (int key : matrix.getRows().keySet()) {

            BigFraction y = new BigFraction(new BigInteger(String.valueOf(len + 1 - key * 2)), (new BigInteger(String.valueOf(len))));

            for (int col : constraints.getMatrix().getRows().get(key)) {
                BigFraction x = new BigFraction(new BigInteger(String.valueOf(len + 1 - col * 2)), (new BigInteger(String.valueOf(len))));
                sequences.put(reverseKaos(x, y, constraints.getWordlength()), new int[]{key, col});
            }


        }


    }

    /**
     *
     * @param s Sting representing a codeword
     * @return row position in mCGR
     */
    public int getRow(String s) {

        return this.sequences.get(s)[0];
    }

    /**
     *
     * @param s Sting representing a codeword
     * @return col position in mCGR
     */
    public int getCol(String s) {

        return this.sequences.get(s)[1];
    }

    public HashMap<String, int[]> getSequences() {
        return this.sequences;
    }


    /**
     * Method to save codewords
     *
     * @param path Output path (as string)
     */
    public void saveAsDNA(String path) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(path));
            for (String sequence : this.getSequences().keySet()) {
                int row = this.getRow(sequence);
                int col = this.getCol(sequence);
                writer.write(">" + "row-" + row + "," + "col-" + col + "\n");
                writer.write(sequence + "\n");
            }
            writer.flush();
            writer.close();
            System.out.println("The translation is done, constrained DNA is saved in " + path);
        } catch (IOException e) {
            System.out.println("Can not find output path, your data wont be stored. ");
        }

    }

    /**
     * @param x      last x-coordinate of a CGR
     * @param y      last y-coordinate of a CGR
     * @param length sequence length
     */
    private static String reverseKaos(BigFraction x, BigFraction y, int length) {
        String sequence;
        BigFraction zx;
        BigFraction zy;
        BigFraction one = new BigFraction(1);
        BigFraction minusOne = new BigFraction(-1);
        StringBuilder stringBuilder = new StringBuilder();

        for (int i = 0; i < length; i++) {


            if (!isNegative(x) && !isNegative(y)) {
                stringBuilder.append("G");
                zx = one;
                zy = one;
            } else if (!isNegative(x) && isNegative(y)) {
                stringBuilder.append("A");
                zx = one;
                zy = minusOne;
            } else if (isNegative(x) && !isNegative(y)) {
                stringBuilder.append("C");
                zx = minusOne;
                zy = one;
            } else {
                stringBuilder.append("T");
                zx = minusOne;
                zy = minusOne;
            }

            x = x.multiply(new BigFraction(2)).subtract(zx);
            y = y.multiply(new BigFraction(2)).subtract(zy);

        }
        sequence = stringBuilder.reverse().toString();
        return sequence;
    }


}
