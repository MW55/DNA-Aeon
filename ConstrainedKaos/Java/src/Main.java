
import constraint.*;
import cgr.ReverseKaos;


/**
 * ConstrainedKaos
 *
 * @author Hannah Franziska Löchel
 */

public class Main {
    private static int length = 0;
    private static int hp = -1;
    private static boolean mode = false;
    static String fileInput = "";
    private static double gc;
    private static double gcStart = 0;
    private static double gcEnd = 1;
    private static int plotSize = 0;
    private static String fileOutput = "";
    private static boolean out = false;
    private static boolean constrains = false;
    private static boolean gcCon = false;
    private static Input input;

    public static void main(String[] args) {

        if (args.length > 0) {
            if (args.length % 2 == 0) {

                for (int i = 0; i < args.length; ) {
                    setParameters(args[i], args[i + 1]);
                    i = i + 2;
                }


                if (length > 0 && out) {


                    Constraints constraintsKaos = calculateConstrains();
                    GCContent gcContent = calculateGC();
                    constraintsKaos.filterGC(gcContent);
                    System.out.println("Lexicographic condition is " + mode + ".");
                    if (!mode) {
                        ReverseKaos reverseKaos = new ReverseKaos(constraintsKaos);
                        reverseKaos.saveAsDNA(fileOutput);
                    } else {
                        Link links = new Link(length, input.getInput());
                        constraintsKaos.filterMotifs(links);
                        ReverseKaos reverseKaos2 = new ReverseKaos(constraintsKaos);
                        reverseKaos2.saveAsDNA(fileOutput);


                    }


                    double ratio = constraintsKaos.ratio();
                    System.out.println("Ratio of allowed sequences: " + ratio * 100 + " %");


                    if (plotSize > 0) {
                        Plot.plot(constraintsKaos.getMatrix(), plotSize);
                    }


                } else {
                    System.out.println("Sequence length and destination of output file are both required.");
                }
            } else {
                System.out.println("Wrong number of input parameters. Try -help for further information.");
            }

        } else help();


    }

    private static Constraints calculateConstrains() {
        if (constrains) {
            input = calculateInput(hp, fileInput);
            if (input.getInput().keySet().size() > 0) {
                return new Constraints(length, input.getInput());
            }

        } else {
            input = new Input();
        }
        return new Constraints(length);

    }

    private static Input calculateInput(int hp, String fileInput) {
        if (hp > 0 && !fileInput.equals("") && hp <= length) {
            return new Input(fileInput, length).add(new Input(hp));
        } else if (hp > 0 && hp <= length) {
            return new Input(hp);
        } else return new Input(fileInput, length);

    }

    private static GCContent calculateGC() {
        if (!gcCon) {
            return new GCContent(length, gcStart, gcEnd);
        } else {
            return new GCContent(length, gc);
        }
    }


    private static void setParameters(String arg, String arg1) {
        switch (arg) {
            //required
            case "-length":
                try {
                    length = Integer.parseInt(arg1);
                } catch (NumberFormatException ex) {
                    System.out.println("Length has to be a positive natural number.");
                    System.exit(1);
                }
                if (length > 12) {
                    System.out.println("Warning: the process might be terminated by the JVM, if not enough RAM is provided.");
                }
                break;
            case "-output":
                fileOutput = arg1;
                out = true;
                break;
            //at least one required 2 for gc interval
            case "-hp":
                try {
                    hp = Integer.parseInt(arg1);
                } catch (NumberFormatException ex) {
                    System.out.println("The number of homopolymers must be a positive natural number.");
                    System.exit(1);
                }
                if (hp > 0) {
                    constrains = true;
                }
                break;
            case "-input":
                fileInput = arg1;
                constrains = true;
                break;
            case "-gc":
                try {
                    gc = Double.parseDouble(arg1);
                } catch (NumberFormatException ex) {
                    System.out.println("GC content has to be a floating point number.");
                    System.exit(1);
                }
                if (gc > 1 | gc < 0) {
                    System.out.println("GC content has to be addressed as floating point number between 0 to 1.");
                    System.exit(1);
                }

                gcCon = true;
                break;
            case "-gcStart":
                try {
                    gcStart = Double.parseDouble(arg1);
                } catch (NumberFormatException ex) {
                    System.out.println("GC content has to be a floating point number.");
                    System.exit(1);
                }
                if (gcStart > 1 | gcStart < 0) {
                    System.out.println("GC content has to be addressed as floating point number between 0 to 1.");
                    System.exit(1);
                }
                break;
            case "-gcEnd":
                try {
                    gcEnd = Double.parseDouble(arg1);
                } catch (NumberFormatException ex) {
                    System.out.println("GC content has to be a floating point number.");
                    System.exit(1);
                }
                if (gcEnd > 1 | gcEnd < 0) {
                    System.out.println("GC content has to be addressed as floating point number between 0 to 1.");
                    System.exit(1);
                }
                break;
            //optional
            case "-plot":
                try {
                    plotSize = Integer.parseInt(arg1);
                } catch (NumberFormatException ex) {
                    System.out.println("Plot size has to be a natural number.");
                    System.exit(1);
                }
                break;
            //optional
            case "-lex":
                try {
                    mode = Boolean.parseBoolean(arg1);
                } catch (IllegalArgumentException ex) {
                    System.out.println("Lex has to be true or false, the program will start without lexicographic condition");

                }
                break;

            case "-help":
                help();
            default:
                break;

        }
    }

    private static void help() {
        System.out.println(
                "ConstrainedKaos Version 1.1.1\n" +
                        "Arguments:\n" +
                        "Required:\n" +
                        "-length: length of desired codewords (necessary)\n" +
                        "-output: destination to save codewords\n" +
                        "At least one argument is required\n" +
                        "-hp: length of the homopolymer, which should be constrained\n" +
                        "-input: path to fasta file with constrained sequences\n" +
                        "-gc:  gc content as float\n" +
                        "or for an interval both of the following are required\n" +
                        "-gcStart: GC content start as float\n" +
                        "-gcEnd: GC content end as float\n" +
                        "optional:\n" +
                        "-plot: size as integer of the dots (we recommend 1 - 5) in the mCGR plot, if -plot is not used, no plot will be created\n" +
                        "-lex: default is false to generate all code words, true generates code words for lexicographic encoding \n"

        );
    }


}