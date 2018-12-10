/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.cnv.preprocess;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.rtg.bed.BedHeader;
import com.rtg.bed.BedReader;
import com.rtg.bed.BedRecord;
import com.rtg.bed.BedWriter;
import com.rtg.util.QuickSort;
import com.rtg.util.QuickSortDoubleIntProxy;
import com.rtg.util.StringUtils;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.intervals.SequenceNameLocusSimple;

/**
 * Holds regions grouped by chromosome
 */
public class RegionDataset {

  /**
   * Reads a dataset from a BED file importing all columns
   * @param file the bed file
   * @return the dataset
   * @throws IOException if there was a problem reading the input
   */
  public static RegionDataset readFromBed(File file) throws IOException {
    try (BedReader br = BedReader.openBedReader(null, file, 0)) {
      final String[] columnNames = getColumnNames(br.getHeader());
      final RegionDataset dataset = new RegionDataset(columnNames);
      loadBedRecords(br, dataset, true /* add new columns */);
      Diagnostic.userLog("Read dataset containing " + dataset.size() + " regions and " + dataset.columns() + " columns from " + file);
      dataset.integrity();
      return dataset;
    }

  }

  /**
   * Reads a dataset from a BED file
   * @param file the bed file
   * @param desiredColumns restrict column loading to the provided columns
   * @return the dataset
   * @throws IOException if there was a problem reading the input
   */
  public static RegionDataset readFromBed(File file, List<Column> desiredColumns) throws IOException {
    try (BedReader br = BedReader.openBedReader(null, file, 0)) {
      final List<String> columnNames = Arrays.asList(getColumnNames(br.getHeader()));
      for (Column c : desiredColumns) {
        if (!columnNames.contains(c.getName())) {
          throw new IOException("Input file: " + file + " does not contain the expected column: " + c.getName());
        }
      }
      final RegionDataset dataset = new RegionDataset(desiredColumns);
      dataset.setColumnIndex(columnNames);

      loadBedRecords(br, dataset, desiredColumns == null);
      Diagnostic.userLog("Read dataset containing " + dataset.size() + " regions and " + dataset.columns() + " columns from " + file);
      dataset.integrity();
      return dataset;
    }
  }

  private static void loadBedRecords(BedReader br, RegionDataset dataset, boolean addNewColumns) throws IOException {
    String sequenceName = null;
    while (br.hasNext()) {
      final BedRecord rec = br.next();
      if (addNewColumns && dataset.size() == 0 && rec.getAnnotations().length > dataset.columns()) {
        dataset.identityIndex(rec.getAnnotations().length);
        for (int i = 0; i < rec.getAnnotations().length; ++i) {
          dataset.addColumn();
        }
      }
      final String newSequence = rec.getSequenceName();
      if (sequenceName == null || !sequenceName.equals(newSequence)) {
        sequenceName = newSequence;
        Diagnostic.userLog("Reading data for sequence " + sequenceName);
      }
      dataset.add(sequenceName, rec.getStart(), rec.getEnd(), rec.getAnnotations());
    }
  }

  // Look for the last header row containing tab-separated values as being column labels
  private static String[] getColumnNames(BedHeader header) {
    String[] names = new String[0];
    final String[] bh = header.getHeaderLines();
    for (int i = bh.length - 1; i >= 0; --i) {
      final String line = bh[i];
      if (line.length() > 0 && line.charAt(0) == BedWriter.COMMENT_CHAR && line.indexOf('\t') != -1) {
        final String[] parts = StringUtils.split(line.substring(1), '\t');
        if (parts.length > 3) {
          names = Arrays.copyOfRange(parts, 3, parts.length);
          break;
        }
      }
    }
    return names;
  }


  private RegionColumn mRegions;
  private final ArrayList<Column> mColumns = new ArrayList<>();

  private int[] mColumnIndexes; // Contains a lookup when calling add with multiple values, when only a subset of columns is required.


  RegionDataset(String[] desiredColumns) {
    this(Arrays.stream(desiredColumns).map(StringColumn::new).collect(Collectors.toList()));
  }

  RegionDataset(List<Column> desiredColumns) {
    mRegions = new RegionColumn("region");
    mColumns.addAll(desiredColumns);
    identityIndex(mColumns.size());
  }

  /**
   * Construct a new dataset with supplied regions but no data columns.
   * @param regions the region column.
   */
  public RegionDataset(RegionColumn regions) {
    mRegions = new RegionColumn("region");
    for (SequenceNameLocus r : regions.values()) {
      mRegions.add(r);
    }
  }

  public RegionDataset copy() {
    return filter(i -> true);
  }


  /**
   * Sets the names of each column being provided in input data when calling add.
   * @param headerColumns the names of each data column supplied during upcoming add calls.
   */
  void setColumnIndex(List<String> headerColumns) {
    final List<String> desiredColumns = getColumnNames();
    final int[] columnIndexes = new int[desiredColumns.size()];
    int i = 0;
    for (final String columnName : desiredColumns) {
      columnIndexes[i++] = headerColumns.indexOf(columnName);
    }
    mColumnIndexes = columnIndexes;
  }

  /** Sets default column index, assuming add will supply one value for each column, in the correct order. */
  private void identityIndex(int length) {
    mColumnIndexes = ArrayUtils.identity(length);
  }

  /**
   * @return the number of rows in the dataset
   */
  public int size() {
    return mRegions.size();
  }

  /**
   * @return the number of columns in the dataset
   */
  public int columns() {
    return mColumns.size();
  }

  List<String> getColumnNames() {
    return mColumns.stream().map(Column::getName).collect(Collectors.toCollection(ArrayList::new));
  }

  /**
   * Get all of the columns in the dataset. It is permitted to modify this list (e.g. removing or adding columns).
   * @return the columns
   */
  public List<Column> getColumns() {
    return mColumns;
  }

  /**
   * @return the column containing the regions
   */
  public RegionColumn regions() {
    return mRegions;
  }

  /**
   * Gets a specific column
   * @param i the column index. A negative index allows indexing from the end of the dataset (e.g. -1 will return the last column).
   * @return the column
   */
  public Column column(int i) {
    return mColumns.get(i < 0 ? mColumns.size() + i : i);
  }

  StringColumn addColumn() {
    return addColumn(new StringColumn("col_" + columns()));
  }

  /**
   * Add the specified column
   * @param column column to add
   * @param <T> type of the column
   * @return the added column
   */
  public <T extends Column> T addColumn(T column) {
    mColumns.add(column);
    mColumnIndexes = null;
    return column;
  }

  /**
   * Add the specified column at the given position
   * @param pos position to add the column at
   * @param column column to add
   * @param <T> type of the column
   * @return the added column
   */
  public <T extends Column> T addColumn(int pos, T column) {
    mColumns.add(pos, column);
    mColumnIndexes = null;
    return column;
  }

  /**
   * Gets the index for the column with the given name
   * @param name the column name
   * @return the index of that column, or -1 if not found
   */
  public int columnId(String name) {
    for (int i = 0; i < mColumns.size(); ++i) {
      if (columnName(i).equals(name)) {
        return i;
      }
    }
    return -1;
  }

  /**
   * Gets the name of a column
   * @param col the column index
   * @return the column name
   */
  public String columnName(int col) {
    return column(col).getName();
  }

  /**
   * Adds a new data row
   * @param sequenceName the sequence name
   * @param start the start position
   * @param end the end position
   * @param value data column values
   */
  public void add(String sequenceName, int start, int end, double... value) {
    mRegions.add(new SequenceNameLocusSimple(sequenceName, start, end));
    for (int i = 0; i < mColumnIndexes.length; ++i) {
      if (mColumnIndexes[i] < value.length) {
        asNumeric(i).add(value[mColumnIndexes[i]]);
      }
    }
  }

  /**
   * Adds a new data row
   * @param sequenceName the sequence name
   * @param start the start position
   * @param end the end position
   * @param value data column values
   */
  public void add(String sequenceName, int start, int end, String... value) {
    mRegions.add(new SequenceNameLocusSimple(sequenceName, start, end));
    for (int i = 0; i < mColumnIndexes.length; ++i) {
      if (mColumnIndexes[i] < value.length) {
        column(i).add(value[mColumnIndexes[i]]);
      }
    }
  }

  /**
   * Adds a new empty data row at a specific position
   * @param row the row index to insert at
   * @param region the region
   */
  public void add(int row, SequenceNameLocus region) {
    mRegions.add(row, region);
    for (int i = 0; i < mColumns.size(); ++i) {
      column(i).add(row, "");
    }
  }

  /**
   * Removes the specified data row
   * @param row the row index
   */
  public void remove(int row) {
    mRegions.remove(row);
    for (int i = 0; i < mColumns.size(); ++i) {
      column(i).remove(row);
    }
  }

  /**
   * Computes the median of a column
   * @param col the column to operate on
   * @return the median
   */
  public double median(int col) {
    return asNumeric(col).median();
  }

  /**
   * Computes the mean of a column
   * @param col the column to operate on
   * @return the mean
   */
  public double mean(int col) {
    return asNumeric(col).mean();
  }

  /**
   * Computes the median of a column, weighted by region length
   * @param col the column to operate on
   * @return the median
   */
  public double weightedMedian(final int col) {
    final NumericColumn column = asNumeric(col);
    final double[] values = column.getValues();
    final int[] weights = new int[size()];
    double tot = 0;
    for (int i = 0; i < values.length; ++i) {
      weights[i] = mRegions.get(i).getLength();
      tot += weights[i];
    }
    final double mid = tot / 2;
    QuickSort.sort(new QuickSortDoubleIntProxy(values, weights));
    tot = 0;
    for (int i = 0; i < values.length; ++i) {
      tot += weights[i];
      if (tot >= mid) {
        return values[i];
      }
    }
    return Double.NaN;
  }

  /**
   * Computes the mean of a column, weighted by region length
   * @param col the column to operate on
   * @return the mean
   */
  public double weightedMean(int col) {
    final NumericColumn column = asNumeric(col);
    double sum = 0;
    double sumLengths = 0;
    for (int i = 0; i < size(); ++i) {
      final int length = mRegions.get(i).getLength();
      sum += column.get(i) * length;
      sumLengths += length;
    }
    return sum / sumLengths;
  }


  void integrity() {
    for (Column c : mColumns) {
      assert c.size() == size();
    }
  }

  /**
   * Gets a column of data, ensuring that it is numeric. If this column has not previously been
   * accessed as numeric data, it is parsed, replacing the original unparsed column. Any
   * unparsable values are treated as Double.NaN
   * @param col the column to get
   * @return the numeric column
   */
  public NumericColumn asNumeric(int col) {
    final Column c = column(col);
    if (c instanceof NumericColumn) {
      return (NumericColumn) c;
    }
    if (c instanceof StringColumn) {
      final NumericColumn c2 = new NumericColumn(c.getName(), "%g");
      c2.set((StringColumn) c);
      mColumns.set(col < 0 ? mColumns.size() + col : col, c2); // Overwrite original
      return c2;
    }
    throw new UnsupportedOperationException("Can't convert column: " + c.getClass() + " to numeric");
  }

  /**
   * Create a dataset only containing rows accepted by the supplied predicate
   * @param p the filtering predicate
   * @return the filtered dataset
   */
  public RegionDataset filter(RegionPredicate p) {
    final RegionDataset d = new RegionDataset(new String[]{});
    d.mRegions = (RegionColumn) mRegions.filter(p);
    for (int col = 0; col < columns(); ++col) {
      d.addColumn(column(col).filter(p));
    }
    return d;
  }

  /**
   * Writes out the dataset in BED format
   * @param bw the BED writer
   * @throws IOException if a problem occurs during writing
   */
  public void write(BedWriter bw) throws IOException {
    bw.writeComment("chrom\tstart\tend\t" + StringUtils.join("\t", getColumnNames()));
    for (int i = 0; i < size(); ++i) {
      bw.write(getBedRecord(i));
    }
    Diagnostic.userLog("Wrote dataset containing " + size() + " regions and " + columns() + " columns");
  }

  BedRecord getBedRecord(int i) {
    final SequenceNameLocus r = regions().get(i);
    final String[] annots = new String[columns()];
    for (int c = 0; c < columns(); ++c) {
      annots[c] = column(c).toString(i);
    }
    return new BedRecord(r.getSequenceName(), r.getStart(), r.getEnd(), annots);
  }
}
