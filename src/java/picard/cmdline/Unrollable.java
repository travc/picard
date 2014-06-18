package picard.cmdline;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Iterables;
import htsjdk.samtools.util.Lazy;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * A little bit of magic!
 * @author mccowan
 */
final public class Unrollable<T> extends ArrayList<String> {
    final Class<T> clazz;
    private Unrollable(final Class<T> clazz) {
        this.clazz = clazz;
    }

    final Lazy<List<T>> unrolled = new Lazy<List<T>>(new Lazy.LazyInitializer<List<T>>() {
        @Override
        public List<T> make() {
            return reduce(Unrollable.this.clazz, Unrollable.this);
        }
    });
    
    /**  */
    public List<T> unrolled() {
        return unrolled.get();
    }

    public static <T> Unrollable<T> of(final Class<T> clazz) {
        return new Unrollable<T>(clazz);
    }
    
    private static <T> List<T> reduce(final Class<T> clazz, final List<String> arguments) {
        final List<String> serializedValues;
        if (!allExtantFilePaths(arguments)) {
            serializedValues = arguments;
        } else {
            serializedValues = new ArrayList<String>();
            for (final String path : arguments) {
                try {
                    serializedValues.addAll(IOUtil.slurpLines(new File(path)));
                } catch (FileNotFoundException e) {
                    throw new RuntimeIOException(e);
                }
            }
        }
        final List<T> values = new ArrayList<T>();
        for (final String serializedValue : serializedValues) {
            values.add((T) CommandLineParser.constructFromString(clazz, serializedValue));    
        }
        return values;
    }
    
    private static boolean allExtantFilePaths(final Collection<String> paths) {
        boolean allExtantFiles = true;
        for (final String argument : paths) {
            allExtantFiles &= argument != null;
            if (!allExtantFiles) break;
            try {
                final File f = new File(argument);
                allExtantFiles = f.isFile();
            } catch (final Exception e) {
                allExtantFiles = false;
            }
        }
        return allExtantFiles;
    }
}
