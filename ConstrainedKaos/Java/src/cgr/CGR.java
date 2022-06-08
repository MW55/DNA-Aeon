package cgr;

import org.apache.commons.math3.fraction.BigFraction;

/**
 * Interface for CGR and reverse CGR
 * ConstrainedKaos  20.11.19
 *
 * @author Hannah Franziska Löchel
 */
public abstract class CGR {

    /**
     * reuturns -, if input is negative
     *
     * @param bigFraction input bigFraction
     * @return true for negative, false for positive
     */
    static boolean isNegative(BigFraction bigFraction) {
        char c = bigFraction.toString().toCharArray()[0];
        return c == '-';
    }

}
