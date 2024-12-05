from sklearn.metrics import confusion_matrix
from sklearn.svm import SVC
import numpy as np
import xgboost as xgb
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
)
from sklearn.preprocessing import LabelEncoder


def percentage_conf_matrix(y_true, y_pred):
    cm = confusion_matrix(y_true, y_pred)
    cm_percentage = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis] * 100
    return cm_percentage


def get_results(train_data, test_data1, test_data2, model):
    X_train = train_data.drop(columns=["curation_status"])
    y_train = train_data["curation_status"]

    X_test1 = test_data1.drop(columns=["curation_status"])
    y_test1 = test_data1["curation_status"]

    X_test2 = test_data2.drop(columns=["curation_status"])
    y_test2 = test_data2["curation_status"]

    # Align test data with training features
    X_test_aligned1 = X_test1[X_train.columns]
    X_test_aligned1 = X_test_aligned1.reindex(columns=X_train.columns, fill_value=0)

    X_test_aligned2 = X_test2[X_train.columns]
    X_test_aligned2 = X_test_aligned2.reindex(columns=X_train.columns, fill_value=0)

    label_encoder = LabelEncoder()
    y_train_encoded = label_encoder.fit_transform(y_train)

    y_test_encoded1 = label_encoder.transform(y_test1)
    y_test_encoded2 = label_encoder.transform(y_test2)

    if isinstance(model, SVC) or isinstance(model, xgb.XGBClassifier):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test_aligned1 = scaler.transform(X_test_aligned1)
        X_test_aligned2 = scaler.transform(X_test_aligned2)

    model = model()

    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    cv_scores = cross_val_score(
        model, X_train, y_train_encoded, cv=cv, scoring="accuracy"
    )

    model.fit(X_train, y_train_encoded)

    y_test_pred1 = model.predict(X_test_aligned1)
    y_test_pred2 = model.predict(X_test_aligned2)

    test_accuracy1 = accuracy_score(y_test_encoded1, y_test_pred1)
    test_precision1 = precision_score(y_test_encoded1, y_test_pred1, average="weighted")
    test_recall1 = recall_score(y_test_encoded1, y_test_pred1, average="weighted")
    test_f1_score1 = f1_score(y_test_encoded1, y_test_pred1, average="weighted")

    if hasattr(model, "predict_proba"):
        auc_roc1 = roc_auc_score(
            y_test_encoded1, model.predict_proba(X_test_aligned1)[:, 1]
        )
    else:
        auc_roc1 = "N/A"

    conf_matrix1 = percentage_conf_matrix(y_test_encoded1, y_test_pred1)

    test_accuracy2 = accuracy_score(y_test_encoded2, y_test_pred2)
    test_precision2 = precision_score(y_test_encoded2, y_test_pred2, average="weighted")
    test_recall2 = recall_score(y_test_encoded2, y_test_pred2, average="weighted")
    test_f1_score2 = f1_score(y_test_encoded2, y_test_pred2, average="weighted")

    if hasattr(model, "predict_proba"):
        auc_roc2 = roc_auc_score(
            y_test_encoded2, model.predict_proba(X_test_aligned2)[:, 1]
        )
    else:
        auc_roc2 = "N/A"

    conf_matrix2 = percentage_conf_matrix(y_test_encoded2, y_test_pred2)

    results = {
        "cv_mean": cv_scores.mean(),
        "Negatives from experiments": {
            "Accuracy": test_accuracy1,
            "Precision": test_precision1,
            "Recall": test_recall1,
            "F1-score": test_f1_score1,
            "AUC-ROC": auc_roc1,
            "Confusion Matrix": conf_matrix1.tolist(),
        },
        "Random Negatives": {
            "cv_mean": cv_scores.mean(),
            "Accuracy": test_accuracy2,
            "Precision": test_precision2,
            "Recall": test_recall2,
            "F1-score": test_f1_score2,
            "AUC-ROC": auc_roc2,
            "Confusion Matrix": conf_matrix2.tolist(),
        },
    }

    return results
