package constraint;

import java.util.HashSet;

/**
 * ConstrainedKaos  27.06.20
 * Class for GC content restrictions
 *
 * @author Hannah Franziska Löchel and Marius Welzel
 */
public class GCContent {
    private HashSet<Integer> gc;
    private int size;


    /**
     * Constructor for the desired GC content restriction
     *
     * @param length     of desired codewords
     * @param proportion of desired gc content
     */
    public GCContent(int length, double proportion) {
        size = (int) Math.pow(2, length - 1);
        int[] contend = calculateGC(length, size);
        this.gc = new HashSet<>();


        int len = size * 2;
        for (int i = 0; i < contend.length; i++) {

            //Bottom GC
            if (length - length * proportion == (double) contend[i] + 1) {
                this.gc.add(size - i);
            }


            //Top AT
            if (length - length * proportion == (double) contend[i]) {
                this.gc.add(len);
            }


            len--;

        }


    }

    /**
     * Constructor for the desired gc content restriction in an interval
     *
     * @param length of desired codewords
     * @param start  of interval
     * @param end    of interval
     */
    public GCContent(int length, double start, double end) {
        size = (int) Math.pow(2, length - 1);
        int[] contend = calculateGC(length, size);

        this.gc = new HashSet<>();

        int len = size * 2;
        for (int i = 0; i < contend.length; i++) {

            //Bottom GC
            if (length - length * start >= contend[i] + 1 && length - length * end <= contend[i] + 1) {
                this.gc.add(size - i);
            }

            //Top AT
            if (length - length * start >= contend[i] && length - length * end <= contend[i]) {
                this.gc.add(len);
            }

            len--;
        }


    }

    /**
     * calculation of gc content
     *
     * @param length of desired codeword
     * @param size   half size of the CGR-matrix
     * @return array with total GC content
     */
    private int[] calculateGC(int length, int size) {


        int[] contend = new int[size];

        for (int i = 0; i < length; i++) {
            int switcher = 1;
            int pos = (int) Math.pow(2, i);
            for (int j = 0; j < size; j++) {

                if (j % pos == 0) {
                    switcher = 1 - switcher;
                }
                contend[j] = contend[j] + switcher;
            }
        }
        return contend;


    }

    /**
     * @return Hashset with columns with desired GC contend in CGR matrix
     */
    public HashSet<Integer> getGc() {
        return this.gc;
    }

    /**
     * @return length of CGR matrix
     */
    public int getLen() {
        return this.size * 2;
    }


}
