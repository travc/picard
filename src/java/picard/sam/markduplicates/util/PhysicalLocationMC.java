package picard.sam.markduplicates.util;

/**
 * @author nhomer
 */
/** Stores the minimal information needed for optical duplicate detection. */
public class PhysicalLocationMC implements OpticalDuplicateFinder.PhysicalLocation {

    // Information used to detect optical dupes
    short readGroup = -1;
    short tile = -1;
    short x = -1, y = -1;
    short libraryId;

    public PhysicalLocationMC(final OpticalDuplicateFinder.PhysicalLocation rec) {
        this.setReadGroup(rec.getReadGroup());
        this.setTile(rec.getTile());
        this.setX(rec.getX());
        this.setY(rec.getY());
        this.setLibraryId(rec.getLibraryId());
    }

    public short getReadGroup() { return this.readGroup; }
    public void  setReadGroup(final short rg) { this.readGroup = rg; }
    public short getTile() { return this.tile; }
    public void  setTile(final short tile) { this.tile = tile; }
    public short getX() { return this.x; }
    public void  setX(final short x) { this.x = x; }
    public short getY() { return this.y; }
    public void  setY(final short y) { this.y = y;}
    public short getLibraryId() { return this.libraryId; }
    public void  setLibraryId(final short libraryId) { this.libraryId = libraryId; }
}