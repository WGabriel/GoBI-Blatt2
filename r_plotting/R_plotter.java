package r_plotting;

import java.io.*;
import java.util.Map;

public class R_plotter implements Runnable {

    public static void main(String[] args) {
        System.out.println("Main started.");
        File file = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt2\\bigOutput\\read.mappinginfo");

        MappingInfo mi = MappingInfo.getMappingInfoFromFile(file);

        // -----Create and execute Barplot R-Command-----
        String rCom = "png(\"C:/Users/Gabriel/Desktop/GoBI/Blatt2/r_plots/barplot.png\", width = 1200, height = 548)\n"
                + "x <- c(" + mi.allReads + ", " + mi.nonSplitReads + ", " + mi.nonSplitReadsWithoutMismatch
                + ", " + mi.splitReads + "," + mi.splitReadsWithoutMismatch + ", " + mi.splitReadsWithoutMismatch5bp + ")\n"
                + "names <- c(\"all reads\",\"non-split reads\",\"non-split reads w/o mismatches\",\"split-reads\","
                + "\"split-reads w/o mismatches\",\"split-reads w/o mismatches\\n and regions >= 5bp\")\n"
                + "cols <- c(\"cadetblue\", \"aquamarine4\",\"aquamarine3\",\"coral4\",\"coral3\",\"coral2\")\n"
                + "barplot(x, names.arg=names, col=cols)\n"
                + "grid(nx=NA, ny=NULL)\n"
                + "par(new=TRUE)\n barplot(x, names.arg=names, col=cols)\n" //re-print barplot OVER the grid
                + "title(main=\"Number of reads in diverse categories\",ylab=\"Number of reads\")\n"
                + "box()\n"
                + "dev.off()\n";
        File rScript = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt2\\r_plots\\myTempScript.R");
        writeStringToFile(rCom, rScript);
        runRscript(rScript);

        // -----Create and execute Mutations Distribution R-Command-----
        // mutations2read: 0->106, 1->178, 2->125, 3->59, 4->23, 5->7, 6->2
        rCom = "png(\"C:/Users/Gabriel/Desktop/GoBI/Blatt2/r_plots/distribution_mutations.png\", width = 1200, height = 548)\n"
                + "plot(1, 1, type = \"n\", ylim=c(0,3000000),yaxp=c(0,3000000,10),xlim=c(0,12),xaxp=c(0,12,12),ann=F,panel.first=grid())\n"
                + "title(main=\"Empirical Distribution Function of mutations per read\",xlab=\"Mutations in read\",ylab=\"Number of reads\")\n";
        for (Map.Entry<Integer, Integer> entry : mi.mutations2reads.entrySet()) {
            rCom += "points(" + entry.getKey() + "," + entry.getValue() + ",col=\"blue\",pch=16)\n";
        }
        rCom += "dev.off()\n";
        rScript = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt2\\r_plots\\myTempScript.R");
        writeStringToFile(rCom, rScript);
        runRscript(rScript);

        // -----Create and execute FragmentLength Distribution R-Command-----
        rCom = "png(\"C:/Users/Gabriel/Desktop/GoBI/Blatt2/r_plots/distribution_fragmentLength.png\", width = 1200, height = 548)\n"
                + "plot(1, 1, type = \"n\", ylim=c(0,50000),yaxp=c(0,50000,20),xlim=c(0,500),xaxp=c(0,500,20),ann=F,panel.first=grid())\n"
                + "title(main=\"Empirical Distribution Function of fragment length per read\",xlab=\"Fragment length in read\",ylab=\"Number of reads\")\n";
        for (Map.Entry<Integer, Integer> entry : mi.fragLength2reads.entrySet()) {
            rCom += "points(" + entry.getKey() + "," + entry.getValue() + ",col=\"blue\",pch=\".\")\n";
        }
        rCom += "dev.off()\n";
        rScript = new File("C:\\Users\\Gabriel\\Desktop\\GoBI\\Blatt2\\r_plots\\myTempScript.R");
        writeStringToFile(rCom, rScript);
        runRscript(rScript);

    }

    private static void runRscript(File script) {
        try {
            System.out.println("R-Script: " + script.getAbsolutePath());
            Process p = new ProcessBuilder("Rscript", "\"" + script.getAbsolutePath() + "\"").start();
            p.waitFor();
            String line;
            // --- stdout ---
            BufferedReader stdout = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((line = stdout.readLine()) != null)
                System.out.println("StdOut: " + line);
            // --- sterr ---
            BufferedReader stderr = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            while ((line = stderr.readLine()) != null)
                System.err.println("StdErr: " + line);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void writeStringToFile(String string, File file) {
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(file));
            writer.write(string);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    @Override
    public void run() {
    }

}

