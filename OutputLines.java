import java.util.ArrayList;

public class OutputLines {
    ArrayList<String> readsFw = new ArrayList<>();
    ArrayList<String> readsRw = new ArrayList<>();
    ArrayList<String> readMappingInfo = new ArrayList<>();

    public OutputLines(ArrayList<String> readsFw, ArrayList<String> readsRw, ArrayList<String> readMappingInfo) {
        this.readsFw = readsFw;
        this.readsRw = readsRw;
        this.readMappingInfo = readMappingInfo;
    }
}
