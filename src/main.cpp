#include "Args/Args.h"
#include "Logger/Logger.h"

int main(int argc, char *argv[]) {
  Logger logger = Logger();
  Args args(logger);
  args.parse_args(argc, argv);
  args.debug_args();
  logger.log(Logger::INFO, "Further logs will be in the file specified by the --o flag.");

  logger = Logger(args.output_file_path);
  logger.info("This will be logged into the file.");

  return 0;
}
