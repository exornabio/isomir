#include <iostream>
#include "config.h"
#include "detection.h"

int main(int argc, char *argv[]) {
    Config *config = new Config();
    config->parse_arg(argc, argv);
    Detection detection = Detection(config);
    detection.detect();
    detection.clear();
}
