RESULT := $(shell mkdir -p $(BUILD_DIR))

all: $(BUILD_DIR)/arepoconfig.h $(BUILD_DIR)/hdf5_util.h

$(BUILD_DIR)/arepoconfig.h: $(CONFIG) $(AREPO_ROOT)/prepare-config.py
	$(PYTHON) $(AREPO_ROOT)/prepare-config.py $(CONFIG) $(BUILD_DIR)

$(BUILD_DIR)/hdf5_util.h:
	cp $(AREPO_ROOT)/src/hdf5_util.h $(BUILD_DIR)
