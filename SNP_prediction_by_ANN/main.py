from data_preprocessing import load_and_preprocess_data
from dataset import DNADataset
from model import DNAClassifier
from train_eval import train_model, evaluate_model, plot_confusion_matrix
import torch
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix


def main():
    # 1. Load and preprocess data
    X, y, feature_names, seq_len = load_and_preprocess_data()
    print(f"Feature vector shape: {X.shape}")
    # 2. Train/test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
    # 3. Dataset and DataLoader
    train_dataset = DNADataset(X_train, y_train)
    test_dataset = DNADataset(X_test, y_test)
    train_loader = DataLoader(train_dataset, batch_size=16, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=16, shuffle=False)
    # 4. Model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = DNAClassifier(input_dim=X.shape[1]).to(device)
    criterion = torch.nn.CrossEntropyLoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    # 5. Train
    print("\nTraining model...")
    train_model(model, train_loader, criterion, optimizer, device, epochs=30)
    # 6. Evaluate
    print("\nEvaluating model...")
    preds, trues = evaluate_model(model, test_loader, device)
    acc = accuracy_score(trues, preds)
    print(f"Test Accuracy: {acc*100:.2f}%")
    print("\nClassification Report:")
    print(classification_report(trues, preds))
    cm = confusion_matrix(trues, preds)
    plot_confusion_matrix(cm, classes=["Class 0", "Class 1"])

if __name__ == "__main__":
    main() 