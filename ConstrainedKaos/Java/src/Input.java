import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

/**
 * ConstrainedKaos   11.11.19
 *
 * @author Hannah Franziska Löchel
 */
public class Input {
    private HashMap<Integer, HashSet<String>> input;


    /**
     * Returns input for hp constraint
     *
     * @param hp length of homopolymer
     */
    public Input(int hp) {
        this.input = new HashMap<>();
        HashSet<String> hpSet = new HashSet<>();
        char[] dna = {'a', 'g', 't', 'c'};
        for (char base : dna) {
            char[] con = new char[hp];
            Arrays.fill(con, base);
            String constrain = new String(con);
            hpSet.add(constrain);
        }
        this.input.put(hp, hpSet);

    }

    public Input() {
        this.input = new HashMap<>();
    }

    /**
     * Combines two inputs with different constarints
     *
     * @param input2 to combine
     * @return combined input
     */
    public Input add(Input input2) {
        for (int length : input2.getInput().keySet()) {
            if (this.input.containsKey(length)) {
                HashSet<String> l = this.input.get(length);
                l.addAll(input2.getInput().get(length));
            } else {
                HashSet<String> l = input2.getInput().get(length);
                this.input.put(length, l);
            }
        }


        return this;
    }

    /**
     * Input for constrained sequences in fasta format
     *
     * @param name path
     */
    public Input(String name, int length) {

        this.input = new HashMap<>();

        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(name));
            String line;
            String key = "";
            String sequence;
            StringBuilder stringBuilder = new StringBuilder();

            while ((line = bufferedReader.readLine()) != null) {
                //
                if (!key.equals("") && Pattern.matches("[AGCTagct]*", line)) {
                    stringBuilder.append(line);
                    sequence = stringBuilder.toString().toLowerCase().trim().replace("\n", "").replace("\r", "");

                    //input.put(key, sequence.toLowerCase());
                    if (sequence.length() <= length) {
                        if (this.input.containsKey(sequence.length())) {
                            this.input.get(sequence.length()).add(sequence);
                        } else {
                            HashSet<String> hashSet = new HashSet<>();
                            hashSet.add(sequence);
                            this.input.put(sequence.length(), hashSet);
                        }

                    } else {
                        System.out.println(">" + key + "\n" + sequence + "\nis longer then codeword length and will be dismissed.");
                    }


                } else if (line.startsWith(">")) {
                    key = line.substring(1);
                    stringBuilder = new StringBuilder();


                } else if (Pattern.matches("[AGCTagct]*", line)) {

                    stringBuilder.append(line);

                }


            }


            bufferedReader.close();

        } catch (IOException e) {
            System.out.println("Could not find input file");
            System.exit(1);
        }


    }

    /**
     * Returns the mCGR of the input with sequence length as key and sequences as values
     *
     * @return HashMap with sequence length as key and sequences as values
     */
    public HashMap<Integer, HashSet<String>> getInput() {
        return this.input;
    }


}
