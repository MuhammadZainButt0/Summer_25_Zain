import logging
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score

 
def model(df):
 try:

    x = df.drop(columns = ["labels","Protein ID"])
    y = df["labels"]

    
    X_train, X_test, Y_train, Y_test = train_test_split(x,y,test_size=0.2,random_state=30)
    logging.info(f"\n\nSelected Training Data:\n {X_train}")
    logging.info(f"\n\nSelected Testing Data:\n {X_test}")
    

    p = Pipeline([
        ("scaler", StandardScaler()),
        # ("SVR",SVR(kernel="linear")),
        ( "Random_Forest", RandomForestClassifier(n_estimators=100, random_state=30))
    ])
    p.fit(X_train, Y_train)


    # scaler = StandardScaler()
    # scaler.fit(X_train)
    # X_train_scaled = scaler.transform(X_train)

    # model = RandomForestClassifier(n_estimators=100, random_state=42)
    # model.fit(X_train_scaled,Y_train)

    Y_pred = p.predict(X_test)
    logging.info(f"\nPrediction: {Y_pred}")
    arry = []
    for i in Y_test:
        arry.append(i)
    logging.info(f"\nActual: {arry}")


    logging.info(f"\nAccuracy: {accuracy_score(Y_test,Y_pred)}")
    logging.info(f"\nConfusion Matrix:\n {confusion_matrix(Y_test,Y_pred)}")
    logging.info(f"\nClassification Report:\n {classification_report(Y_test,Y_pred)}")
    

    return
 except:
    print(f"Model Implementation Failed")











    
    

