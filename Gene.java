import java.util.HashMap;

public class Gene {

    public String gene_id;
    // String = transcript Ids
    public HashMap<String, Transcript> transcripts = new HashMap<String, Transcript>();

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
        this.gene_id = new String(another.gene_id);
    }


}
