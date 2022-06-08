package constraint;

import cgr.Kaos;
import org.apache.commons.math3.fraction.BigFraction;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * ConstrainedKaos  19.06.20
 * Class to construct and resize mCGR
 *
 * @author Hannah Franziska Löchel
 */

public class mCGR {
    /**
     * mCGR in HashMap structure
     */
    private HashMap<Integer, HashSet<Integer>> rows;
    /**
     * length of mCGR
     */
    private int len;

    mCGR() {
        this.rows = new HashMap<>();
        this.len = 0;
    }

    mCGR(HashSet<String> constrains, int zoom) {
        this.len = (int) Math.pow(2, zoom);

        this.rows = new HashMap<>();

        for (String constrain : constrains) {
            Kaos kaos = new Kaos(constrain);
            int col = toMatrix(kaos.getxLast(), this.len);
            int row = toMatrix(kaos.getyLast(), this.len);
            if (this.rows.containsKey(row)) {
                HashSet<Integer> cols = this.rows.get(row);
                cols.add(col);
                this.rows.put(row, cols);
            } else {
                HashSet<Integer> cols = new HashSet<>();
                cols.add(col);
                this.rows.put(row, cols);
            }
        }

    }

    /**
     * Constructor for mCGR with one constraint
     *
     * @param constrain as String
     */
    mCGR(String constrain) {
        Kaos kaos = new Kaos(constrain);
        this.len = (int) Math.pow(2, constrain.length());
        int col = toMatrix(kaos.getxLast(), this.len);
        int row = toMatrix(kaos.getyLast(), this.len);
        rows = new HashMap<>();
        HashSet<Integer> cols = new HashSet<>();
        cols.add(col);
        rows.put(row, cols);
    }

    /**
     * Function to convert BigInteger fraction CGR to matrix coordinates
     *
     * @param last BigInteger Fraction
     * @param len  CGR Matrix length
     * @return matrix coordinate
     */
    private int toMatrix(BigFraction last, int len) {
        return (int) Math.ceil((last.doubleValue() + 1) * (len >> 1));
    }

    /**
     * Constructor for mCGR of constrained homopolymers
     *
     * @param hp length of constrained homopolymers
     */
    mCGR(int hp) {
        char[] con = new char[hp];
        Arrays.fill(con, 'a');
        String constrain = new String(con);
        Kaos kaos = new Kaos(constrain);
        this.len = (int) Math.pow(2, hp);
        int col = toMatrix(kaos.getxLast(), this.len);
        int row = toMatrix(kaos.getyLast(), this.len);
        int end = this.len;
        rows = new HashMap<>();
        HashSet<Integer> cols = new HashSet<>();
        cols.add(row);
        cols.add(col);
        rows.put(col, cols);
        rows.put(end, cols);
    }

    /**
     * Constructor for mCGR of GC content, allowed sequences are stored, a fully allowed row is marked with -1
     *
     * @param gc GCContend object
     */
    mCGR(GCContent gc) {
        HashSet<Integer> init = new HashSet<>();
        rows = new HashMap<>();
        init.add(-1);
        for (int element : gc.getGc()) {
            rows.put(element, init);
        }
        this.len = gc.getLen();
    }


    /**
     * Constructor for mCGR
     *
     * @param matrix as HashMap of mCGR with codewords or constraints
     * @param length size of mCGR
     */
    mCGR(HashMap<Integer, HashSet<Integer>> matrix, int length) {
        this.rows = matrix;
        this.len = length;
    }

    mCGR doubleSize() {
        HashMap<Integer, HashSet<Integer>> d = new HashMap<>();


        for (int row : this.rows.keySet()) {
            int newRow = row * 2;
            int newRowNext = row * 2 - 1;
            HashSet<Integer> newColOne;
            newColOne = new HashSet<>();
            HashSet<Integer> newColTwo;
            newColTwo = new HashSet<>();
            HashSet<Integer> cols = this.rows.get(row);

            for (int col : cols) {
                newColOne.add(col * 2);
                newColOne.add(col * 2 - 1);
                newColTwo.add(col * 2);
                newColTwo.add(col * 2 - 1);
            }
            d.put(newRow, newColOne);
            d.put(newRowNext, newColTwo);
        }
        return new mCGR(d, this.len * 2);
    }

    /**
     * Method to tile over the next iteration of CGR
     *
     * @param mCGR for tiling
     * @return changed mCGR
     */
    mCGR tiling(mCGR mCGR) {

        for (int row : mCGR.rows.keySet()) {
            int newRow = row + mCGR.len;
            HashSet<Integer> newColOne;

            if (rows.containsKey(row)) {
                newColOne = this.rows.get(row);
            } else {
                newColOne = new HashSet<>();
                this.rows.put(row, newColOne);
            }

            HashSet<Integer> newColTwo;

            if (this.rows.containsKey(newRow)) {
                newColTwo = this.rows.get(newRow);
            } else {
                newColTwo = new HashSet<>();
                this.rows.put(newRow, newColTwo);
            }

            HashSet<Integer> cols = mCGR.rows.get(row);

            for (int col : cols) {
                newColOne.add(col);
                newColOne.add(col + mCGR.len);
                newColTwo.add(col);
                newColTwo.add(col + mCGR.len);
            }
        }


        return this;
    }


    /**
     * Returns the mCGR as HashMap, where the keys are the rows (y-values) and the values are the columns stored in a HashSet (y-Values)
     *
     * @return mCGR as HashMap/sparsematrix
     */
    public HashMap<Integer, HashSet<Integer>> getRows() {
        return this.rows;
    }

    mCGR add(mCGR mCGR) {
        for (int row : mCGR.rows.keySet()) {
            if (this.rows.containsKey(row)) {
                HashSet<Integer> col = this.rows.get(row);
                col.addAll(mCGR.rows.get(row));
            } else {
                this.rows.put(row, mCGR.rows.get(row));
            }
        }

        return this;
    }

    /**
     * Filters the allowed sequences of a mCGR containing the forbidden sequences with a given gcContend
     *
     * @param gcContent to filter with
     * @return filterd mCGR sparsematrix with codewords
     */
    mCGR filter(GCContent gcContent) {
        HashMap<Integer, HashSet<Integer>> rows2 = new HashMap<>();

        for (int row : gcContent.getGc()) {
            HashSet<Integer> init = new HashSet<>();

            for (int i = 0; i < this.len; i++) {
                init.add(i + 1);
            }

            if (this.rows.containsKey(row)) {
                init.removeAll(this.rows.get(row));
            }

            rows2.put(row, init);
        }
        this.rows = rows2;
        return this;
    }

    /**
     * @return length of mCGR sparsematrix
     */
    public int getLen() {
        return this.len;
    }

    void setLen(int len) {
        this.len = (int) Math.pow(2, len);
    }

    void spotStart(mCGR preC) {
        int iterator = this.getLen() / preC.getLen();
        for (int row : preC.rows.keySet()) {
            for (int i = 1; i < iterator; i++) {
                int newRow = row + preC.getLen() * i;
                HashSet<Integer> newColOne;

                if (rows.containsKey(row)) {
                    newColOne = this.rows.get(row);
                } else {
                    newColOne = new HashSet<>();
                    this.rows.put(row, newColOne);
                }

                HashSet<Integer> cols = preC.rows.get(row);

                for (int col : cols) {
                    newColOne.add(col);
                    newColOne.add(col + preC.getLen() * i);
                }
                this.rows.put(newRow, newColOne);
            }
        }
    }

    public mCGR filter(Link link) {
        for (mCGR m: link.getMCGRs()) {
            for (int row: m.getRows().keySet()) {
                HashSet<Integer> col =  m.getRows().get(row);
                for (int column: col    ) {
                    if(this.getRows().containsKey(row)){
                        this.getRows().get(row).remove(column);
                    }
                }
            }
        }
        return this;
    }
}
