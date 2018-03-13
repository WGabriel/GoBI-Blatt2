import java.util.HashMap;

public class Gene {
    public String gene_id;
    // Key = transcript_Ids
    public HashMap<String, Transcript> transcripts = new HashMap<>();

    // Default Constructor
    public Gene() {
        super();
    }

    // Constructor
    public Gene(String gene_id, HashMap<String, Transcript> transcripts) {
        super();
        this.gene_id = gene_id;
        this.transcripts = transcripts;
    }

    // Constructor for cloning a Gene
    public Gene(Gene another) {
        //Deep copy of Exon list required
        for (HashMap.Entry<String, Transcript> t : another.transcripts.entrySet()) {
            this.transcripts.put(t.getKey(), new Transcript(t.getValue()));
        }
        this.gene_id = another.gene_id;
    }


}
