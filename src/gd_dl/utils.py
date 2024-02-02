import argparse


def set_prep(
		mode,
		valid_idx,
		train_val_ratio
	):
	train_s = list()
	val_s = list()

	if mode == 'valid':
		val_s = [l.strip() for l in open('./data/valid.set')]
		return train_s, val_s

	elif mode == 'train':
		i = 0
		with open('./data/train.set') as f:
			for line in f:
				if i % train_val_ratio != valid_idx:
					train_s.append(line.strip())
				else:
					val_s.append(line.strip())
				i+=1
		return train_s, val_s


def str2bool(v):
	if v.lower() in ['yes', 'true', 't', 'y', '1']:
		return True
	elif v.lower() in ['no', 'false', 'f', 'n', '0']:
		return False
	else:
		raise argparse.ArgumentTypeError('Boolean value expected')
