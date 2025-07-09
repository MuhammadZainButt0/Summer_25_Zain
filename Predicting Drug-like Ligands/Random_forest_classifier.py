import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score,precision_recall_curve,mean_squared_error


def main():
    features = pd.read_csv("Descriptors_encoding.csv")


    x = features.drop(["Drug", "compound"], axis=1)
    y = features["Drug"]

    X_train, X_test, Y_train, Y_test = train_test_split(x,y,test_size=0.2,random_state=40)


    model = RandomForestClassifier(n_estimators=100, random_state=40)
    model.fit(X_train,Y_train)
    Y_pred = model.predict(X_test)
    print("Predicted values:",Y_pred)
    arry = []
    for i in Y_test:
        arry.append(i)
    print("Actual Value:",arry)

    print("Confusion Matrix:\n",confusion_matrix(Y_test,Y_pred))
    print("\nClassification Report:\n",classification_report(Y_test,Y_pred))
    print("\n Accuracy score:",accuracy_score(Y_test,Y_pred))
    print("\n Precision score:",precision_recall_curve(Y_test,Y_pred))
    print("\n Mean Square Error:", mean_squared_error(Y_test, Y_pred))



if __name__ == "__main__":
    main()

