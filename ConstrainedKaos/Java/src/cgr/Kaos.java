package cgr;

import org.apache.commons.math3.fraction.BigFraction;


/**
 * ConstrainedKaos  01.11.19
 * Class to construct CGR
 *
 * @author Hannah Franziska Löchel
 */
public class Kaos extends CGR {


    private BigFraction[] c = {new BigFraction(1), new BigFraction(-1)};
    private BigFraction[] g = {new BigFraction(-1), new BigFraction(-1)};
    private BigFraction[] a = {new BigFraction(-1), new BigFraction(1)};
    private BigFraction[] t = {new BigFraction(1), new BigFraction(1)};


    private BigFraction xLast;
    private BigFraction yLast;


    /**
     * CGR on sequence
     *
     * @param dna as String
     */
    public Kaos(String dna) {

        char[] data = dna.toCharArray();
        BigFraction sf = new BigFraction(1, 2);
        BigFraction x = new BigFraction(0);
        BigFraction y = new BigFraction(0);


        for (char i : data) {
            x = x.add(sf.multiply(getCoord(i)[0].subtract(x)));
            y = y.add(sf.multiply(getCoord(i)[1].subtract(y)));
        }

        this.xLast = x;
        this.yLast = y;

    }

    public BigFraction getxLast() {
        return this.xLast;
    }

    public BigFraction getyLast() {
        return this.yLast;
    }


    /**
     * returns coordinates of CGR
     *
     * @param element dna base
     * @return coordinates
     */
    private BigFraction[] getCoord(char element) {
        if (element == 'a') {
            return a;
        } else if (element == 'g') {
            return g;
        } else if (element == 'c') {
            return c;
        }
        return t;
    }


}
