<?php
// submit_registration.php

if (isset($_POST['register'])) {
    // Retrieve form data
    $year = $_POST['year'];
    $month = $_POST['month'];
    $name_course = $_POST['name_course'];
    $name_participant = $_POST['name_participant'];
    $email_participant = $_POST['email_participant'];
    $name_PI = $_POST['name_PI'];
    $email_PI = $_POST['email_PI'];
    $project_number = $_POST['project_number'];
    $fund = $_POST['fund'];
    $dept_id = $_POST['dept_id'];

    try {
        // Connect to SQLite database
        $db = new PDO('sqlite:courses_participants.db');
        $db->setAttribute(PDO::ATTR_ERRMODE, PDO::ERRMODE_EXCEPTION);

        // Create table if it doesn't exist
        $db->exec("CREATE TABLE IF NOT EXISTS participants (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            year INTEGER,
            month TEXT,
            name_course TEXT,
            name_participant TEXT,
            email_participant TEXT,
            name_PI TEXT,
            email_PI TEXT,
            project_number TEXT,
            fund TEXT,
            dept_id TEXT
        )");

        // Prepare SQL statement to prevent SQL injection
        $stmt = $db->prepare("INSERT INTO participants (
            year, month, name_course, name_participant, email_participant,
            name_PI, email_PI, project_number, fund, dept_id
        ) VALUES (
            :year, :month, :name_course, :name_participant, :email_participant,
            :name_PI, :email_PI, :project_number, :fund, :dept_id
        )");

        // Bind parameters
        $stmt->bindParam(':year', $year, PDO::PARAM_INT);
        $stmt->bindParam(':month', $month, PDO::PARAM_STR);
        $stmt->bindParam(':name_course', $name_course, PDO::PARAM_STR);
        $stmt->bindParam(':name_participant', $name_participant, PDO::PARAM_STR);
        $stmt->bindParam(':email_participant', $email_participant, PDO::PARAM_STR);
        $stmt->bindParam(':name_PI', $name_PI, PDO::PARAM_STR);
        $stmt->bindParam(':email_PI', $email_PI, PDO::PARAM_STR);
        $stmt->bindParam(':project_number', $project_number, PDO::PARAM_STR);
        $stmt->bindParam(':fund', $fund, PDO::PARAM_STR);
        $stmt->bindParam(':dept_id', $dept_id, PDO::PARAM_STR);

        // Execute statement
        $stmt->execute();

        // Close the database connection
        $db = null;

        // Redirect to thank you page
        header("Location: thank_you.html");
        exit();

    } catch (PDOException $e) {
        echo "An error occurred while saving your registration: " . $e->getMessage();
    }
} else {
    // If the form wasn't submitted properly, redirect back to the form
    header("Location: registration.html");
    exit();
}
?>