/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.vcf.validator;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;

import com.rtg.launcher.AbstractCli;
import com.rtg.util.Resources;
import com.rtg.util.StringUtils;
import com.rtg.util.TextTable;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.header.FilterField;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.MetaType;
import com.rtg.vcf.header.VcfHeader;
import com.rtg.vcf.header.VcfNumber;
import com.rtg.vcf.header.VcfNumberType;
import com.rtg.vcf.validator.NumericRuleSet.DoubleConverter;
import com.rtg.vcf.validator.NumericRuleSet.LongConverter;
import com.rtg.vcf.validator.RuleSet.FieldType;
import com.rtg.vcf.validator.RuleSet.StringConverter;

/**
 * A small utility for validating the contents of
 * a VCF file conform to expected value ranges.
 */
public class VcfValidatorCli extends AbstractCli {

  private static final String XRULES = "Xrules";
  private static final String XVERBOSE = "Xverbose";

  private static final int MAX_RECORD_WARNINGS = 10;

  private final Map<String, RuleSet<?>> mInfoRules = new HashMap<>();
  private final Map<String, RuleSet<?>> mFormatRules = new HashMap<>();

  /*
   * Cannot detect this as our reading code interferes with it by removing the PASS
   * filter every time you add a new one. The only case we could possibly detect
   * is if the PASS filter is the last one in the list as the reader does not
   * deal with this.
   */
  //private long mNumberPassFilterWithOthers = 0;

  private long mNumberFieldsNotInHeader = 0;
  private long mNumberInvalidRecords = 0;
  private long mNumberFieldsNotInRules = 0;
  private long mNumberInvalidQualValues = 0;
  private boolean mPrintedRecord = false;
  private boolean mCountedRecord = false;
  private boolean mVerbose;

  @Override
  public String moduleName() {
    return "vcfvalidator";
  }

  @Override
  public String description() {
    return null;
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Validates the contents of a VCF file conform to expected value ranges.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerExtendedHelp();
    mFlags.registerRequired(File.class, "FILE", "VCF format file to be validated").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(XRULES, File.class, "FILE", "File defining rules for validation of VCF input").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(XVERBOSE, "Set to output all failed records to error output instead of the first " + MAX_RECORD_WARNINGS).setCategory(CommonFlagCategories.REPORTING);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File vcfFile = (File) mFlags.getAnonymousValue(0);
    readRules((File) mFlags.getValue(XRULES));
    mVerbose = mFlags.isSet(XVERBOSE);
    int mismatchedHeaders = 0;
    try (final VcfReader reader = VcfReader.openVcfReader(vcfFile)) {
      final VcfHeader header = reader.getHeader();
      final Map<String, FilterField> filters = new HashMap<>();
      for (FilterField field : header.getFilterLines()) {
        filters.put(field.getId(), field);
      }
      final Map<String, InfoField> infos = new HashMap<>();
      for (InfoField field : header.getInfoLines()) {
        infos.put(field.getId(), field);
        final RuleSet<?> rule = mInfoRules.get(field.getId());
        if (rule == null) {
          err.println("VCF INFO field " + field.getId() + " not contained in the specified rule set.");
          ++mNumberFieldsNotInRules;
        } else {
          if (rule.getMetaType() != field.getType()) {
            err.println("VCF INFO field " + field.getId() + " does not match the expected meta-data type.");
            ++mismatchedHeaders;
          }
          if (rule.getVcfNumber().getNumberType() != field.getNumber().getNumberType() || rule.getVcfNumber().getNumber() != field.getNumber().getNumber()) {
            err.println("VCF INFO field " + field.getId() + " does not match the expected number of values.");
            ++mismatchedHeaders;
          }
        }
      }
      final Map<String, FormatField> formats = new HashMap<>();
      for (FormatField field : header.getFormatLines()) {
        formats.put(field.getId(), field);
        final RuleSet<?> rule = mFormatRules.get(field.getId());
        if (rule == null) {
          err.println("VCF FORMAT field " + field.getId() + " not contained in the specified rule set.");
          ++mNumberFieldsNotInRules;
        } else {
          if (rule.getMetaType() != field.getType()) {
            err.println("VCF FORMAT field " + field.getId() + " does not match the expected meta-data type.");
            ++mismatchedHeaders;
          }
          if (rule.getVcfNumber().getNumberType() != field.getNumber().getNumberType() || rule.getVcfNumber().getNumber() != field.getNumber().getNumber()) {
            err.println("VCF FORMAT field " + field.getId() + " does not match the expected number of values.");
            ++mismatchedHeaders;
          }
        }
      }
      if (mismatchedHeaders > 0) {
        final TextTable warnings = new TextTable();
        addWarningToTable(warnings, mismatchedHeaders, "Number of fields in headers that conflict with rule definitions:");
        addWarningToTable(warnings, mNumberFieldsNotInRules, "Number of fields in headers not defined in rules:");
        err.println();
        err.println(warnings.toString());
        return 1;
      }

      while (reader.hasNext()) {
        final VcfRecord current = reader.next();
        if (!(current.getFilters().size() == 1 && VcfRecord.MISSING.equals(current.getFilters().get(0)))) {
          for (String filter : current.getFilters()) {
//            if ("PASS".equals(filter)) {
//              if (current.getFilters().size() > 1) {
//                warning(err, current, "VCF FILTER field PASS present with other FILTER fields.");
//                mNumberPassFilterWithOthers++;
//              }
            if (!"PASS".equals(filter) && !filters.containsKey(filter)) {
               warning(err, current, "VCF FILTER field + " + filter + " not present in the header.");
               ++mNumberFieldsNotInHeader;
            }
          }
        }
        for (String info : current.getInfo().keySet()) {
          if (!infos.containsKey(info)) {
            warning(err, current, "VCF INFO field " + info + " not present in the header.");
            ++mNumberFieldsNotInHeader;
          }
          final RuleSet<?> rule = mInfoRules.get(info);
          if (rule != null) {
            try {
              rule.validateRecord(current);
            } catch (RuleValidationException e) {
              warning(err, current, e.getMessage());
            }
          }
        }
        for (String format : current.getFormats()) {
          if (!formats.containsKey(format)) {
            warning(err, current, "VCF FORMAT field " + format + " not present in the header.");
            ++mNumberFieldsNotInHeader;
          }
          final RuleSet<?> rule = mFormatRules.get(format);
          if (rule != null) {
            try {
              rule.validateRecord(current);
            } catch (RuleValidationException e) {
              warning(err, current, e.getMessage());
            }
          }
        }
        if (!VcfRecord.MISSING.equals(current.getQuality())) {
          Double val;
          try {
            val = Double.parseDouble(current.getQuality());
          } catch (NumberFormatException e) {
            val = null;
          }
          if (val == null || Double.isInfinite(val) || Double.isNaN(val)) {
            warning(err, current, "QUAL value is invalid.");
            ++mNumberInvalidQualValues;
          }
        }
        mPrintedRecord = false;
        mCountedRecord = false;
      }
    }
    final TextTable warnings = new TextTable();

    addWarningToTable(warnings, mNumberFieldsNotInRules, "Number of fields in headers not defined in rules:");
    addWarningToTable(warnings, mNumberFieldsNotInHeader, "Number of fields in records not defined in header:");
    //addWarningToTable(warnings, mNumberPassFilterWithOthers, "Number of records with conflicting PASS filter:");
    addWarningToTable(warnings, mNumberInvalidQualValues, "Number of records with invalid QUAL fields:");
    addWarningToTable(warnings, mNumberInvalidRecords, "Total number of invalid records:");

    if (warnings.numRows() > 0) {
      err.println();
      err.println(warnings.toString());
    }

    return (mNumberFieldsNotInRules > 0 || mNumberInvalidRecords > 0) ? 1 : 0;
  }

  private void addWarningToTable(TextTable table, long count, String title) {
    if (count > 0) {
      table.addRow(title, Long.toString(count));
    }
  }

  private void warning(PrintStream err, VcfRecord rec, String message) {
    if (!mPrintedRecord) {
      if (!mCountedRecord) {
        ++mNumberInvalidRecords;
        mCountedRecord = true;
      }
      if (!mVerbose && mNumberInvalidRecords > MAX_RECORD_WARNINGS) {
        return;
      }
      err.println();
      err.println(rec.toString());
      mPrintedRecord = true;
    }
    err.println(message);
  }

  private void readRules(File rulesFile) throws IOException {
    final Properties pr = new Properties();
    try (InputStream rulesStream = getRulesStream(rulesFile)) {
      pr.load(rulesStream);
    }

    final Map<String, Map<String, String>> infoDefinitions = new HashMap<>();
    final Map<String, Map<String, String>> formatDefinitions = new HashMap<>();
    for (String property : pr.stringPropertyNames()) {
      final String[] parts = StringUtils.split(property, '.');
      if (parts.length != 3) {
        throw new IOException("Invalid property name: " + property);
      }
      Map<String, String> properties;
      if ("INFO".equals(parts[0])) {
        properties = infoDefinitions.get(parts[1]);
        if (properties == null) {
          properties = new HashMap<>();
          infoDefinitions.put(parts[1], properties);
        }
      } else if ("FORMAT".equals(parts[0])) {
        properties = formatDefinitions.get(parts[1]);
        if (properties == null) {
          properties = new HashMap<>();
          formatDefinitions.put(parts[1], properties);
        }
      } else {
        throw new IOException("Invalid field type: " + parts[0]);
      }
      properties.put(parts[2], pr.getProperty(property));
    }
    for (Entry<String, Map<String, String>> info : infoDefinitions.entrySet()) {
      mInfoRules.put(info.getKey(), createRuleSet(FieldType.INFO, info.getKey(), info.getValue()));
    }
    for (Entry<String, Map<String, String>> format : formatDefinitions.entrySet()) {
      mFormatRules.put(format.getKey(), createRuleSet(FieldType.FORMAT, format.getKey(), format.getValue()));
    }
  }

  private InputStream getRulesStream(File rulesFile) throws IOException {
    final InputStream rulesStream;
    if (rulesFile == null) {
      rulesStream = Resources.getResourceAsStream("com/rtg/vcf/validator/rules.properties");
      if (rulesStream == null) {
        throw new IOException("Unable to load validation rules resource.");
      }
    } else {
      rulesStream = new FileInputStream(rulesFile);
    }
    return rulesStream;
  }

  private RuleSet<?> createRuleSet(FieldType type, String name, Map<String, String> rules) throws IOException {
    final String fieldType = rules.get("TYPE");
    if (fieldType == null) {
      throw new IOException("Rules for " + type + " field " + name + " missing TYPE value.");
    }
    final MetaType metaType = MetaType.parseValue(fieldType);
    if (metaType == null) {
      throw new IOException("Rules for " + type + " field " + name + " has unrecognised TYPE value " + fieldType);
    }
    final String fieldNum = rules.get("NUM");
    if (fieldNum == null) {
      throw new IOException("Rules for " + type + " field " + name + " missing NUM value.");
    }
    final VcfNumber number;
    try {
      number = new VcfNumber(fieldNum);
    } catch (NumberFormatException e) {
      throw new IOException("Rules for " + type + " field " + name + " has invalid NUM field value.");
    }
    if (number.getNumber() < 0 && number.getNumberType() == VcfNumberType.INTEGER) {
      throw new IOException("Rules for " + type + " field " + name + " has negative value for NUM field value.");
    }
    try {
      final RuleSet<?> returnRules;
      switch (metaType) {
        case INTEGER:
        case FLOAT:
          returnRules = createNumericRuleSet(name, type, number, metaType, rules);
          break;
        case FLAG:
          if (type == FieldType.FORMAT) {
            throw new IOException("Rules for FORMAT field " + name + " are set to TYPE value FLAG.");
          }
          if (number.getNumber() != 0) {
            throw new IOException("Rules for INFO field " + name + " has TYPE FLAG and NUM not equal to 0.");
          }
          return new RuleSet<>(name, type, number, metaType, null);
        default:
          returnRules = new RuleSet<>(name, type, number, metaType, new StringConverter());
          break;
      }
      if (rules.containsKey("ENUM")) {
        returnRules.addEnumerationRule(rules.get("ENUM"));
      }
      return returnRules;
    } catch (RuleValidationException e) {
      throw new IOException("Rules for " + type + " field " + name + " has one or more invalid values.");
    }
  }

  private NumericRuleSet<?> createNumericRuleSet(String name, FieldType type, VcfNumber number, MetaType metaType, Map<String, String> rules) throws RuleValidationException {
    final NumericRuleSet<?> returnRules;
    if (metaType == MetaType.FLOAT) {
      returnRules = getFloatRuleSet(name, type, number, metaType);
    } else {
      returnRules = getLongRuleSet(name, type, number, metaType);
      if (!(rules.containsKey("ALLOWINF") && Boolean.parseBoolean(rules.get("ALLOWINF")))) {
        returnRules.addInfinityRule();
      }
      if (!(rules.containsKey("ALLOWNAN") && Boolean.parseBoolean(rules.get("ALLOWNAN")))) {
        returnRules.addNaNRule();
      }
    }
    if (rules.containsKey("GTE")) {
      returnRules.addGreaterThanEqualRule(rules.get("GTE"));
    }
    if (rules.containsKey("GT")) {
      returnRules.addGreaterThanRule(rules.get("GT"));
    }
    if (rules.containsKey("LTE")) {
      returnRules.addLessThanEqualRule(rules.get("LTE"));
    }
    if (rules.containsKey("LT")) {
      returnRules.addLessThanRule(rules.get("LT"));
    }
    return returnRules;
  }

  private NumericRuleSet<Long> getLongRuleSet(String name, FieldType type, VcfNumber number, MetaType metaType) {
    return new NumericRuleSet<>(name, type, number, metaType, new LongConverter());
  }

  private NumericRuleSet<Double> getFloatRuleSet(String name, FieldType type, VcfNumber number, MetaType metaType) {
    return new NumericRuleSet<>(name, type, number, metaType, new DoubleConverter());
  }

  /**
   * Main Running Method
   * @param args the command line arguments
   */
  public static void main(String[] args) {
    System.exit(new VcfValidatorCli().mainInit(args, System.out, System.err));
  }

}
