import os
import torch
import torch.nn.functional as F
from torch_geometric.loader import DataLoader
from torch.utils.data.dataloader import default_collate
from torch_geometric.nn import GATConv, Linear, global_mean_pool
import pandas as pd
from torch import nn
import numpy as np
import random
from datetime import datetime
from asd_final_dataset import (
    GACorih_v3_bi_trainsp0_st12_nca_sc_deg_Dataset,
    GACorih_v3_bi_valsp0_st12_nca_sc_deg_Dataset,
    GACorih_v3_bi_testsp0_st12_nca_sc_deg_Dataset
)
from plot_graphs_final import plot_metrics_bi
import matplotlib.pyplot as plt
from torcheval.metrics.functional import binary_auroc, binary_auprc

# =========================
# Ray / Tune 제거됨 (고정 HP)
# =========================
ver_name = "gma3orihv3st12_nca_sc_deg_bi_vw"

start_time = datetime.now()
seed = 5155
np.random.seed(seed=seed)
random.seed(a=seed)
torch.manual_seed(seed)

# 고정 하이퍼파라미터 (요청대로)
max_epoch   = int(1000)
hid_dim     = int(12)
b_size      = int(1)
layer_count = int(3)
num_heads   = int(4)
fixed_lr    = float(0.0006)

num_gene = int(1)
precision_degree = [10, 20, 50, 100]
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print('device', device)

train_data = GACorih_v3_bi_trainsp0_st12_nca_sc_deg_Dataset(root='./data/bapdata/')
print(train_data[0].x.shape, train_data[0].y.shape)
print(len(train_data))
print('train_data', train_data)
# 메모리 관리
# del train_data  # (아래서 다시 쓰므로 삭제하지 않음)

dev_data = GACorih_v3_bi_valsp0_st12_nca_sc_deg_Dataset(root='./data/bapdata/')
print('dev_data', dev_data)
# del dev_data  # (아래서 다시 쓰므로 삭제하지 않음)


class GCN(torch.nn.Module):
    def __init__(self, hidden_channels, layer_count, num_heads):
        super().__init__()
        self.conv1 = GATConv(1, hidden_channels, heads=num_heads)
        self.conv2 = nn.ModuleList([
            GATConv(in_channels=hidden_channels * num_heads,
                    out_channels=hidden_channels * num_heads)
            for i in range(layer_count)
        ])
        if hidden_channels == 1:
            hidden_layer_size = 1
        else:
            hidden_layer_size = int(hidden_channels * num_heads / 2)
        self.lin1 = Linear(hidden_channels * num_heads, hidden_layer_size, weight_initializer='glorot')
        self.lin2 = Linear(hidden_layer_size, int(hidden_layer_size / 2), weight_initializer='glorot')
        self.lin3 = Linear(int(hidden_layer_size / 2), 1, weight_initializer='glorot')

    def forward(self, data):  # , batch
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.gelu(x)
        attn_weights = None
        for l in self.conv2:
            x, attn_weights = l(x, edge_index, return_attention_weights=True)
            x = F.gelu(x)
        # readout
        x = global_mean_pool(x, data.batch)
        x = F.gelu(self.lin1(x))
        x = F.gelu(self.lin2(x))
        out = self.lin3(x)
        out = torch.sigmoid(out)
        return out, attn_weights


class EarlyStopping:
    def __init__(self, patience=200):
        self.loss = np.inf
        self.patience = 0
        self.patience_limit = patience

    def step(self, loss):
        if self.loss > loss:
            self.loss = loss
            self.patience = 0
        else:
            self.patience += 1

    def is_stop(self):
        return self.patience >= self.patience_limit


def gma3orihv3st12_nca_sc_deg_bi_vw():
    # === 고정 HP/데이터로 순차 학습 ===
    early_stop = EarlyStopping(patience=200)
    #
    train_loader = DataLoader(
        train_data, batch_size=b_size, shuffle=True,
        collate_fn=lambda x: default_collate(x).to(device)
    )
    dev_loader = DataLoader(
        dev_data, batch_size=b_size, shuffle=True,
        collate_fn=lambda x: default_collate(x).to(device)
    )
    #
    model = GCN(hid_dim, layer_count, num_heads).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=fixed_lr)
    criterion = nn.BCELoss().to('cpu')
    train_loss_list = []
    dev_loss_list = []
    acc_list_train = []
    acc_list_dev = []
    acc_list_test = []
    auc_list_train = []
    auc_list_dev = []
    auc_list_test = []
    auprc_list_train = []
    auprc_list_dev = []

    best_train_acc = float(0)
    best_dev_acc = float("-inf")
    best_train_auc = float(0)
    best_dev_auc = float("-inf")

    # 저장 경로(원본 tune.checkpoint_dir 대체)
    save_dir = "."
    ckpt_path = os.path.join(save_dir, "checkpoint.pt")  # (model_state, optim_state) tuple 저장

    for epoch in range(max_epoch):
        model.train()
        epoch_loss = 0
        lb_np = torch.empty((0, 1), dtype=torch.int64)
        predict_np = torch.empty((0, 1), dtype=torch.int64)

        for i, data in enumerate(train_loader):
            hypothesis, attn_weights = model(data.to(device))
            loss = criterion(hypothesis, data.y)

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
            predict_np = torch.cat([predict_np, hypothesis.cpu()], dim=0)
            lb_np = torch.cat([lb_np, data.y.cpu()], dim=0)

        predicted = (predict_np > 0.5).float()
        acc = (predicted == lb_np).float().mean()
        accuracy = acc.item()
        train_auc = binary_auroc(predict_np.reshape(-1), lb_np.reshape(-1)).mean().item()
        auc_list_train.append(train_auc)
        train_loss_val = (epoch_loss / (i + 1))
        train_loss_list.append(train_loss_val)
        acc_list_train.append(accuracy)
        train_auprc = binary_auprc(predict_np.reshape(-1), lb_np.reshape(-1)).mean().item()
        auprc_list_train.append(train_auprc)

        # ====== DEV ======
        model.eval()
        epoch_loss_dev = 0
        lb_np = torch.empty((0, 1), dtype=torch.int64)
        predict_np = torch.empty((0, 1), dtype=torch.int64)
        with torch.no_grad():
            for i, data in enumerate(dev_loader):
                hypothesis, attn_weights = model(data.to(device))
                loss = criterion(hypothesis, data.y)
                epoch_loss_dev += loss.item()
                predict_np = torch.cat([predict_np, hypothesis.cpu()], dim=0)
                lb_np = torch.cat([lb_np, data.y.cpu()], dim=0)

            predicted = (predict_np > 0.5).float()
            dev_acc = (predicted == lb_np).float().mean()
            dev_accuracy = dev_acc.item()
            dev_auc = binary_auroc(predict_np.reshape(-1), lb_np.reshape(-1)).mean().item()
            auc_list_dev.append(dev_auc)
            dev_auprc = binary_auprc(predict_np.reshape(-1), lb_np.reshape(-1)).mean().item()
            auprc_list_dev.append(dev_auprc)

            dev_loss_val = (epoch_loss_dev / (i + 1))
            dev_loss_list.append(dev_loss_val)
            early_stop.step(dev_loss_val)

            if best_dev_acc < dev_accuracy:
                best_dev_acc = dev_accuracy
            if best_train_acc < accuracy:
                best_train_acc = accuracy
            if best_dev_auc < dev_auc:
                best_dev_auc = dev_auc
            if best_train_auc < train_auc:
                best_train_auc = train_auc

            # 원본 tune.checkpoint_dir 대체: 매 epoch 임시 체크포인트 경로
            torch.save((model.state_dict(), optimizer.state_dict()), ckpt_path)

            # 마지막 epoch일 때 결과 저장 로직 유지
            if epoch == (max_epoch - 1):
                # path3 계산 대신 현재 디렉토리 사용
                path3 = "."
                torch.save(model.state_dict(), os.path.join(path3, "model.pt"))
                try:
                    result_df = pd.DataFrame({
                        "loss": train_loss_list, "dev_loss": dev_loss_list,
                        "acc": acc_list_train, "dev_acc": acc_list_dev,
                        "auc": auc_list_train, "dev_auc": auc_list_dev,
                        "auprc": auprc_list_train, "dev_auprc": auprc_list_dev,
                    })
                    result_df.index = np.arange(1, len(result_df) + 1)
                    result_df.to_csv(os.path.join(path3, "result_df.csv"))
                    plt.clf()
                    plot_metrics_bi(result_df)
                    plt.savefig(os.path.join(path3, "baseline_history.png"))
                    plt.close()
                except Exception as e:
                    print("error saving last-epoch artifacts:", e, flush=True)

        if early_stop.is_stop():
            path3 = "."
            torch.save((model.state_dict(), optimizer.state_dict()), ckpt_path)
            torch.save(model.state_dict(), os.path.join(path3, "model.pt"))
            try:
                result_df = pd.DataFrame({
                    "loss": train_loss_list, "dev_loss": dev_loss_list,
                    "acc": acc_list_train, "dev_acc": acc_list_dev,
                    "auc": auc_list_train, "dev_auc": auc_list_dev,
                    "auprc": auprc_list_train, "dev_auprc": auprc_list_dev,
                })
                result_df.index = np.arange(1, len(result_df) + 1)
                result_df.to_csv(os.path.join(path3, "result_df.csv"))
                plt.clf()
                plot_metrics_bi(result_df)
                plt.savefig(os.path.join(path3, "baseline_history.png"))
                plt.close()
            except Exception as e:
                print("error saving early-stop artifacts:", e, flush=True)
            break

    # 학습 종료 후 반환값은 원본에서 tune가 처리하던 것을 대체할 필요 없음
    return ckpt_path


# ====== 실행 (학습) ======
ckpt_path = gma3orihv3st12_nca_sc_deg_bi_vw()

# ====== TEST ======
test_data = GACorih_v3_bi_testsp0_st12_nca_sc_deg_Dataset(root='./data/bapdata/')
test_loader = DataLoader(test_data, batch_size=b_size,
                         collate_fn=lambda x: default_collate(x).to(device))

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
best_trained_model = GCN(hid_dim, layer_count, num_heads).to(device)
optimizer = torch.optim.Adam(best_trained_model.parameters(), lr=fixed_lr)
criterion = nn.BCELoss()

# checkpoint 로드 (원본의 best_trial.checkpoint 대체)
if os.path.exists(ckpt_path):
    model_state, optimizer_state = torch.load(ckpt_path, map_location=device)
    best_trained_model.load_state_dict(model_state)
    try:
        optimizer.load_state_dict(optimizer_state)
    except Exception:
        pass

num_gene = 1
batch_size = b_size
epoch_loss = 0
lb_np = torch.empty((0, 1), dtype=torch.int64)
predict_np = torch.empty((0, 1), dtype=torch.int64)
best_trained_model.eval()
with torch.no_grad():
    for i, data in enumerate(test_loader):
        hypothesis, attn_weights = best_trained_model(data.to(device))
        loss = criterion(hypothesis, data.y)
        epoch_loss += loss.item()
        predict_np = torch.cat([predict_np, hypothesis.cpu()], dim=0)
        lb_np = torch.cat([lb_np, data.y.cpu()], dim=0)
    predicted = (predict_np > 0.5).float()
    test_acc = (predicted == lb_np).float().mean()
    test_accuracy = test_acc.item()
    test_auc = binary_auroc(predict_np.reshape(-1), lb_np.reshape(-1)).mean().item()
    test_auprc = binary_auprc(predict_np.reshape(-1), lb_np.reshape(-1)).mean().item()

print("Best trial test set loss: {:.4f}".format(epoch_loss / (i + 1)))
print("Best trial test set acc: {:.4f}".format(test_accuracy))
print("Best trial test set auc: {:.4f}".format(test_auc))
print("Best trial test set auprc: {:.4f}".format(test_auprc))

torch.save(best_trained_model.state_dict(), ver_name + '.pt')

end_time = datetime.now()
print(end_time)
print(end_time - start_time)







