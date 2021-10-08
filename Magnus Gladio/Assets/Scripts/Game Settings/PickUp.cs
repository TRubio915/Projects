using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PickUp : MonoBehaviour
{
    public enum ItemType
    {
        HealthPotion,
        StaminaPotion
    }

    private Inventory inventory;
    public GameObject itemButton;
    public ItemType itemType;

    // Start is called before the first frame update
    void Start()
    {
        inventory = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<Inventory>();
    }

    private void OnTriggerEnter2D(Collider2D collision)
    {
        if(collision.CompareTag("Player"))
        {
            for(int i = 0; i < inventory.slots.Length; i++)
            {
                if(inventory.isFull[i] == false)
                {
                    // Pick up item
                    Instantiate(itemButton, inventory.slots[i].transform, false);
                    Destroy(gameObject);
                    inventory.isFull[i] = true;                    
                    break;

                }
            }
        }
    }

    public void DisableAnimator()
    {
        GetComponent<Animator>().enabled = false;
    }
}
