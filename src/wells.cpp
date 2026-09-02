#pragma once

struct Well {
	enum Type {Injector, Producer};
	int cell;
	Type type;
	double p_bhp;
};