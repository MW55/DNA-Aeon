package constraint;


import java.util.HashMap;
import java.util.HashSet;


/**
 * ConstrainedKaos  24.03.21
 * Class for mCGRs of undesired links
 *
 * @author Hannah Franziska Löchel
 */
public class Link {
    private HashSet<mCGR> links;


    /**
     * @param length of code words
     * @param input  Input, containing constrained motifs
     */
    public Link(int length, HashMap<Integer, HashSet<String>> input) {
        // key: prepending: value: appending
        this.links = new HashSet<>();
        HashSet<String> prependingCGR = new HashSet<>();
        HashSet<String> appendingCGR = new HashSet<>();


        for (Integer i : input.keySet()) {
            for (String s : input.get(i)) {


                for (int j = 1; j < s.length(); j++) {
                    String prepending = s.substring(0, j);
                    String appending = s.substring(j);
                    //dismiss subsequences combinations longer then codeword length
                    if (prepending.length() <= length && appending.length() <= length) {
                        if (prepending.length() >= appending.length()) {
                            prependingCGR.add(prepending);

                        } else {
                            appendingCGR.add(appending);
                        }
                    }
                }
            }


        }

        for (String pre : prependingCGR) {

            mCGR preC = new mCGR(pre);

            //spot endings
            for (int i = pre.length(); i < length; i++) {
                preC = preC.doubleSize();
            }

            this.links.add(preC);

            for (String ap : appendingCGR) {
                mCGR empty = new mCGR(ap);
                mCGR apC = new mCGR();
                apC.setLen(length);
                apC.spotStart(empty);

                links.add(apC);

            }
        }


    }

    /**
     * returns matrix positions of sequences forming undesired motifs, when concatenated
     *
     * @return HashSet with mCGRs
     */
    public HashSet<mCGR> getMCGRs() {
        return this.links;
    }


}
